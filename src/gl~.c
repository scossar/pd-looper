#include "m_pd.h"
#include <stdlib.h>
#include <math.h>

typedef enum {
  STATE_IDLE,
  STATE_RECORDING,
  STATE_PLAYING
} t_gl_state;

typedef struct _grain {
  int position;
  int ms;
  int samples; // should probably be t_float (depends on sample rate)
} t_grain;

typedef struct _gl {
  t_object x_obj;

  t_grain *x_grains;
  int x_grain_ms;

  int x_input_buffer_ms;
  int x_input_buffer_samples; // t_float instead?
  t_sample *x_input_buffer;

  int x_window_buffer_samples; // t_float?
  t_sample *x_window_buffer;

  int x_write_phase;

  t_gl_state x_state;

  t_float x_sms; // samples per ms
  int x_num_grains;
  t_float x_grain_spread; // some float value to set distance between grains;

  t_float x_mix; // wet/dry

  t_inlet *x_inlet_pos;
  t_float x_f; // dummy arg for MAINSIGNALIN
} t_gl;

static t_class *gl_class = NULL;

static void gl_free(t_gl *x);
static int initialize_buffers(t_gl *x);
static int initialize_grains(t_gl *x);

static void *gl_new(t_floatarg grain_ms, t_floatarg num_grains)
{
  t_gl *x = (t_gl *)pd_new(gl_class);

  x->x_input_buffer_ms = 4000;
  x->x_grain_ms = (grain_ms > 10) ? grain_ms : 10; // todo: find a better
  // default
  x->x_num_grains = (num_grains > 0) ? num_grains : 1;
  if (!initialize_grains(x)) {
    gl_free(x); // I _think_ just Pd will handle the call on return NULL, but...
    return NULL;
  }

  x->x_grain_spread = 0.3f; // testing

  x->x_sms = 0;
  x->x_write_phase = 0;

  x->x_mix = 0.5f;
  x->x_state = STATE_IDLE;

  // initialize with small power of 2 values
  x->x_input_buffer_samples = 1024;
  x->x_window_buffer_samples = 1024;
  if (!initialize_buffers(x)) {
    gl_free(x);
    return NULL;
  }

  x->x_inlet_pos = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);

  outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
}

// TODO: rename, or maybe bave separate methods for initializeing and updating
// grains?
static int initialize_grains(t_gl *x)
{
  x->x_grains = (t_grain *)getbytes(sizeof(t_grain) * x->x_num_grains);
  if (x->x_grains == NULL) {
    pd_error(x, "gl~: failed to allocate memory for grains");
    return 0;
  }

  for (int i = 0; i < x->x_num_grains; i++) {
    t_grain *grain = &x->x_grains[i];
    grain->ms = x->x_grain_ms;
    grain->samples = grain->ms * x->x_sms;
    // rate
    grain->position = 0;
  }

  return 1;
}

int initialize_buffers(t_gl *x)
{
  x->x_input_buffer = getbytes(x->x_input_buffer_samples * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "gl~: unable to allocate memory to input buffer");
    return 0;
  }

  x->x_window_buffer = getbytes(x->x_window_buffer_samples * sizeof(t_sample));
  if (x->x_window_buffer == NULL) {
    pd_error(x, "gl~: unable to allocate memory to window buffer");
    return 0;
  }
  return 1;
}

static void update_input_buffer(t_gl *x)
{
  int buffer_size = 1;
  while (buffer_size < (x->x_input_buffer_ms * x->x_sms)) {
    buffer_size *= 2;
  }

  x->x_input_buffer = (t_sample *)resizebytes(x->x_input_buffer,
                                              x->x_input_buffer_samples * sizeof(t_sample),
                                              buffer_size * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "gl~: unable to resize input buffer");
    return;
  }

  x->x_input_buffer_samples = buffer_size;
  x->x_write_phase = 0;
  post("gl~: (debug) x_input_buffer_samples: %d", x->x_input_buffer_samples);
}

static void update_window_buffer(t_gl *x)
{
  int buffer_size = 1;
  while (buffer_size < (x->x_grain_ms * x->x_sms)) {
    buffer_size *= 2;
  }

  x->x_window_buffer = (t_sample *)resizebytes(x->x_window_buffer,
                                               x->x_window_buffer_samples * sizeof(t_sample),
                                               buffer_size * sizeof(t_sample));

  if (x->x_window_buffer == NULL) {
    pd_error(x, "gl~: unable to resize window buffer");
    return;
  }
  x->x_window_buffer_samples = buffer_size;
  post("gl~: (debug) x_window_buffer_samples: %d", x->x_window_buffer_samples);
}

static void hanning_window(t_gl *x)
{
  int grain_samples = x->x_sms * x->x_grain_ms;
  for (int i = 0; i < grain_samples; i++) {
    float phase = (float)i / (grain_samples - 1);
    x->x_window_buffer[i] = 0.5f * (1.0f - cos(2.0f * M_PI * phase));
  }

  for (int i = 0; i < x->x_num_grains; i++) {
    x->x_grains[i].samples = grain_samples;
    x->x_grains[i].position = 0;
  }

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

static void system_params(t_gl *x, t_float sr)
{
  x->x_sms = sr * 0.001f;
}

static t_int *gl_perform(t_int *w)
{
  t_gl *x = (t_gl *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *in2 = (t_sample *)(w[3]);
  t_sample *out = (t_sample *)(w[4]);
  int n = (int)(w[5]);

  t_sample *input_buffer = x->x_input_buffer;
  int input_buffer_samples = x->x_input_buffer_samples;
  int input_buffer_mask = input_buffer_samples - 1;
  int write_phase = x->x_write_phase;

  t_float g_scale = 1.0f / x->x_num_grains;
  t_float grain_offset = x->x_grain_ms * x->x_sms * x->x_grain_spread;

  while (n--) {
    t_sample f = *in1++;

    if (x->x_state == STATE_RECORDING) {
      input_buffer[write_phase] = f;
    }

    t_sample grain_start = *in2++;

    if (grain_start > 1.0f) grain_start = 1.0f;
    if (grain_start < -1.0f) grain_start = -1.0f;
    grain_start = grain_start * 0.5f + 0.5f; // maybe handle scaling in the
    // patch?
    t_sample grain_output = 0.0f;
    for (int i = 0; i < x->x_num_grains; i++) {
      float full_index = (grain_start * input_buffer_samples) + x->x_grains[i].position + i * grain_offset;
      int grain_index = (int)full_index;
      float frac = full_index - (t_sample)grain_index;
      grain_index = grain_index & input_buffer_mask;
      t_sample grain_sample = cubic_interpolate(input_buffer, grain_index, input_buffer_mask, frac);
      grain_output += grain_sample * g_scale * x->x_window_buffer[x->x_grains[i].position];
      x->x_grains[i].position += 1;
      if (x->x_grains[i].position >= x->x_grains[i].samples) {
        x->x_grains[i].position = 0;
      }
    }

    *out++ = (grain_output * x->x_mix) + (f * (1.0f - x->x_mix));

    write_phase = (write_phase +1) & input_buffer_mask;
  }

  x->x_write_phase = write_phase;
  return (w+6);
}

static void gl_dsp(t_gl *x, t_signal **sp)
{
  dsp_add(gl_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_length);
  system_params(x, sp[0]->s_sr);
  update_input_buffer(x);
  update_window_buffer(x);
  hanning_window(x);
  initialize_grains(x);
}

static void gl_free(t_gl *x)
{
  if (x->x_grains != NULL) {
    freebytes(x->x_grains, x->x_num_grains * sizeof(t_grain));
    x->x_grains = NULL;
  }

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

static void looper_record(t_gl *x)
{
  x->x_state = STATE_RECORDING;
}

static void looper_play(t_gl *x)
{
  x->x_state = STATE_PLAYING;
}

static void resize_window(t_gl *x, t_floatarg f)
{
  x->x_grain_ms = (f > 10) ? f : 10;
  update_window_buffer(x);
  hanning_window(x);
  initialize_grains(x);
}

static void mix(t_gl *x, t_floatarg f)
{
  if (f < 0.0f) f = 0.0f;
  if (f > 1.0f) f = 1.0f;
  x->x_mix = f;
}

static void spread(t_gl *x, t_floatarg f) {
  if (f < 0) f = 0;
  x->x_grain_spread = f;
}

void gl_tilde_setup(void)
{
  gl_class = class_new(gensym("gl~"),
                            (t_newmethod)gl_new,
                            (t_method)gl_free,
                            sizeof(t_gl),
                            CLASS_DEFAULT,
                            A_DEFFLOAT, A_DEFFLOAT, 0);

  class_addmethod(gl_class, (t_method)gl_dsp, gensym("dsp"), A_CANT, 0);
  class_addmethod(gl_class, (t_method)looper_record, gensym("record"), 0);
  class_addmethod(gl_class, (t_method)looper_play, gensym("play"), 0);
  class_addmethod(gl_class, (t_method)resize_window, gensym("resize"), A_FLOAT, 0);
  class_addmethod(gl_class, (t_method)mix, gensym("mix"), A_FLOAT, 0);
  class_addmethod(gl_class, (t_method)spread, gensym("spread"), A_FLOAT, 0);
  CLASS_MAINSIGNALIN(gl_class, t_gl, x_f);
}
