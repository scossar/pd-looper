#include "m_pd.h"
#include <stdlib.h>
#include <math.h>

typedef enum {
  STATE_IDLE,
  STATE_RECORDING,
  STATE_PLAYING
} t_glooper_state;

typedef struct _glooper {
  t_object x_obj;

  int x_input_buffer_ms;
  int x_input_buffer_samples;
  t_sample *x_input_buffer;

  int x_grain_ms;
  int x_grain_samples;
  int x_window_buffer_samples;
  t_sample *x_window_buffer;

  int x_write_phase;
  int x_grain_pos;

  t_glooper_state x_state;

  t_float x_sms; // samples per ms
  int x_num_grains;

  t_float x_mix; // wet/dry

  t_inlet *x_inlet_pos;
  t_float x_f; // dummy arg for MAINSIGNALIN
} t_glooper;

static t_class *glooper_class = NULL;
static int initialize_buffers(t_glooper *x);

static void *glooper_new(t_floatarg grain_ms, t_floatarg num_grains)
{
  t_glooper *x = (t_glooper *)pd_new(glooper_class);

  x->x_input_buffer_ms = 4000;
  x->x_grain_ms = (grain_ms > 10) ? grain_ms : 10; // todo: find a better
  // default
  x->x_num_grains = (num_grains > 0) ? num_grains : 1;

  x->x_sms = 0;
  x->x_write_phase = 0;
  x->x_grain_pos = 0;
  x->x_grain_samples = 0;

  x->x_mix = 0.5f;
  x->x_state = STATE_IDLE;

  // initialize with small power of 2 values
  x->x_input_buffer_samples = 1024;
  x->x_window_buffer_samples = 1024;
  if (!initialize_buffers(x)) {
    return NULL;
  }

  x->x_inlet_pos = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);

  outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
}

int initialize_buffers(t_glooper *x)
{
  x->x_input_buffer = getbytes(x->x_input_buffer_samples * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "glooper~: unable to allocate memory to input buffer");
    return 0;
  }

  x->x_window_buffer = getbytes(x->x_window_buffer_samples * sizeof(t_sample));
  if (x->x_window_buffer == NULL) {
    pd_error(x, "glooper~: unable to allocate memory to window buffer");
    return 0;
  }
  return 1;
}

static void update_input_buffer(t_glooper *x)
{
  int buffer_size = 1;
  while (buffer_size < (x->x_input_buffer_ms * x->x_sms)) {
    buffer_size *= 2;
  }

  x->x_input_buffer = (t_sample *)resizebytes(x->x_input_buffer,
                                              x->x_input_buffer_samples * sizeof(t_sample),
                                              buffer_size * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "glooper~: unable to resize input buffer");
    return;
  }

  x->x_input_buffer_samples = buffer_size;
  x->x_write_phase = 0;
  post("glooper~: (debug) x_input_buffer_samples: %d", x->x_input_buffer_samples);
}

static void update_window_buffer(t_glooper *x)
{
  int buffer_size = 1;
  while (buffer_size < (x->x_grain_ms * x->x_sms)) {
    buffer_size *= 2;
  }

  x->x_window_buffer = (t_sample *)resizebytes(x->x_window_buffer,
                                               x->x_window_buffer_samples * sizeof(t_sample),
                                               buffer_size * sizeof(t_sample));

  if (x->x_window_buffer == NULL) {
    pd_error(x, "glooper~: unable to resize window buffer");
    return;
  }
  x->x_window_buffer_samples = buffer_size;
  post("glooper~: (debug) x_window_buffer_samples: %d", x->x_window_buffer_samples);
}

static void hanning_window(t_glooper *x)
{
  int grain_samples = x->x_sms * x->x_grain_ms;
  for (int i = 0; i < grain_samples; i++) {
    float phase = (float)i / (grain_samples - 1);
    x->x_window_buffer[i] = 0.5f * (1.0f - cos(2.0f * M_PI * phase));
  }
  x->x_grain_samples = grain_samples;
  x->x_grain_pos = 0;
}

// NOTE: for the `inline` directive to be respected by the compiler, probably
// more this method to a header file
// this method is using a look-behind approach to interpolation, try look
// forward, or look both forward and back?
static inline t_sample cubic_interpolate(t_sample *buffer, int phase, int mask, t_sample frac)
{
  t_sample a = buffer[phase];
  t_sample b = buffer[(phase - 1) & mask];
  t_sample c = buffer[(phase - 2) & mask];
  t_sample d = buffer[(phase - 3) & mask];
  t_sample cminusb = c - b;

  return b + frac * (
      cminusb - 0.1666667f * (1.0f - frac) * (
          (d - a - 3.0f * cminusb) * frac + (d + 2.0f * a - 3.0f * b)
      )
  );
}

static void system_params(t_glooper *x, t_float sr)
{
  x->x_sms = sr * 0.001f;
}

static t_int *glooper_perform(t_int *w)
{
  t_glooper *x = (t_glooper *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *in2 = (t_sample *)(w[3]);
  t_sample *out = (t_sample *)(w[4]);
  int n = (int)(w[5]);

  t_sample *input_buffer = x->x_input_buffer;
  int input_buffer_samples = x->x_input_buffer_samples;
  int input_buffer_mask = input_buffer_samples - 1;
  int write_phase = x->x_write_phase;
  int grain_pos = x->x_grain_pos;

  while (n--) {
    t_sample f = *in1++;

    if (x->x_state == STATE_RECORDING) {
      input_buffer[write_phase] = f;
    }

    t_sample gs = *in2++;

    // the expected input for gs will be in the range (-1, 1)
    // scale it to (0, 1)
    if (gs > 1) gs = 1.0f;
    if (gs < -1) gs = -1.0f;
    gs = gs * 0.5f + 0.5f;

    float full_index = (gs * input_buffer_samples) + grain_pos;
    int grain_index = (int)full_index;
    float frac = full_index - grain_index;
    grain_index = grain_index & input_buffer_mask;
    t_sample grain_sample = cubic_interpolate(input_buffer, grain_index, input_buffer_mask, frac);

    grain_sample = input_buffer[grain_index] * x->x_window_buffer[grain_pos];

    *out++ = (grain_sample * x->x_mix) + (f * (1.0f - x->x_mix));

    write_phase = (write_phase +1) & input_buffer_mask;
    grain_pos = grain_pos + 1;
    if (grain_pos >= x->x_grain_samples) grain_pos = 0;
  }

  x->x_grain_pos = grain_pos;
  x->x_write_phase = write_phase;
  return (w+6);
}

static void glooper_dsp(t_glooper *x, t_signal **sp)
{
  dsp_add(glooper_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_length);
  system_params(x, sp[0]->s_sr);
  update_input_buffer(x);
  update_window_buffer(x);
  hanning_window(x);
}

static void glooper_free(t_glooper *x)
{
  if (x->x_input_buffer != NULL) {
    freebytes(x->x_input_buffer, x->x_input_buffer_samples * sizeof(t_sample));
    x->x_input_buffer = NULL;
  }

  if (x->x_window_buffer != NULL) {
    freebytes(x->x_window_buffer, x->x_window_buffer_samples * sizeof(t_sample));
    x->x_window_buffer = NULL;
  }

  if (x->x_inlet_pos != NULL) {
    inlet_free(x->x_inlet_pos);
  }
}

static void looper_record(t_glooper *x)
{
  x->x_state = STATE_RECORDING;
}

static void looper_play(t_glooper *x)
{
  x->x_state = STATE_PLAYING;
}

static void resize_window(t_glooper *x, t_floatarg f)
{
  x->x_grain_ms = (f > 10) ? f : 10;
  update_window_buffer(x);
  hanning_window(x);
}

static void mix(t_glooper *x, t_floatarg f)
{
  if (f < 0.0f) f = 0.0f;
  if (f > 1.0f) f = 1.0f;
  x->x_mix = f;
}

void glooper_tilde_setup(void)
{
  glooper_class = class_new(gensym("glooper~"),
                            (t_newmethod)glooper_new,
                            (t_method)glooper_free,
                            sizeof(t_glooper),
                            CLASS_DEFAULT,
                            A_DEFFLOAT, A_DEFFLOAT, 0);

  class_addmethod(glooper_class, (t_method)glooper_dsp, gensym("dsp"), A_CANT, 0);
  class_addmethod(glooper_class, (t_method)looper_record, gensym("record"), 0);
  class_addmethod(glooper_class, (t_method)looper_play, gensym("play"), 0);
  class_addmethod(glooper_class, (t_method)resize_window, gensym("resize"), A_FLOAT, 0);
  class_addmethod(glooper_class, (t_method)mix, gensym("mix"), A_FLOAT, 0);
  CLASS_MAINSIGNALIN(glooper_class, t_glooper, x_f);
}
