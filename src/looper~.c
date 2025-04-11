// this is turning into a granular synth; not there yet

#include "m_pd.h"
#include <math.h>

typedef enum {
  STATE_IDLE,
  STATE_RECORDING,
  STATE_PLAYING
} t_looper_state;

typedef struct _looper {
  t_object x_obj;

  t_float x_s_per_msec;
  int x_pd_block_size; // possibly not needed

  t_sample *x_input_buffer;
  t_sample *x_window_buffer;
  t_float x_input_buffer_ms; // could be an int?
  int x_input_buffer_samples;

  t_looper_state x_state;

  int x_write_phase; // input buffer write position
  int x_read_phase; // input buffer read position
  int x_loop_pos; // current position in the loop

  int x_loop_start;
  int x_loop_length;
  int x_fade_samples;

  t_float x_f; // dummy arg for CLASS_MAINSIGNALIN
} t_looper;

static t_class *looper_class = NULL;

static void input_buffer_update(t_looper *x);

static void *looper_new(t_floatarg f) {
  t_looper *x = (t_looper *)pd_new(looper_class);

  // used to set the minimum size of x_input_buffer once x->x_x_per_ms is known
  x->x_input_buffer_ms = (f > 1) ? f : 4000.0f;
  x->x_s_per_msec = 0.0f;
  x->x_pd_block_size = 0;

  x->x_state = STATE_IDLE;

  x->x_write_phase = 0;
  x->x_read_phase = 0;
  x->x_loop_start = 0;
  x->x_loop_length = 0;
  x->x_loop_pos = 0;

  x->x_fade_samples = 10 * 64; // hmmm

  x->x_input_buffer_samples = 1024; // initialize to a small power of 2 value
  x->x_input_buffer = getbytes(x->x_input_buffer_samples * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "looper~: unable to assign memory to input buffer");
    return NULL;
  }

  x->x_window_buffer = getbytes(x->x_input_buffer_samples * sizeof(t_sample));
  // handle allocation error here

  x->x_window_buffer = getbytes(x->x_input_buffer_samples * sizeof(t_sample));
  // if (x->x_window_buffer == NULL) {
  //   pd_error(x, "looper~: unable to assign memory to window buffer");
  //   return NULL;
  // }

  outlet_new(&x->x_obj, &s_signal);

  return (void *)x;
}

static void input_buffer_update(t_looper *x)
{
  int buffer_size = 1;
  // note: in the delay implementations I've been (unnecessarily?) adding x_pd_block size here
  while (buffer_size < (x->x_input_buffer_ms * x->x_s_per_msec)) {
    buffer_size *= 2;
  }

  x->x_input_buffer = (t_sample *)resizebytes(x->x_input_buffer,
                                              x->x_input_buffer_samples * sizeof(t_sample),
                                              buffer_size * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "looper~: unable to resize input buffer");
    return;
  }
  x->x_window_buffer = (t_sample *)resizebytes(x->x_input_buffer,
                                              x->x_input_buffer_samples * sizeof(t_sample),
                                              buffer_size * sizeof(t_sample));
  if (x->x_input_buffer == NULL) {
    pd_error(x, "looper~: unable to resize input buffer");
    return;
  }

  x->x_input_buffer_samples = buffer_size;
  x->x_write_phase = 0;
  post("looper~: (debug) x_input_buffer_samples: %d", x->x_input_buffer_samples);
}

static void set_system_params(t_looper *x, int blocksize, t_float sr)
{
  x->x_pd_block_size = blocksize;
  x->x_s_per_msec = sr * 0.001f;
}

static t_int *looper_perform(t_int *w)
{
  t_looper *x = (t_looper *)(w[1]);
  t_sample *in1 = (t_sample *)(w[2]);
  t_sample *out = (t_sample *)(w[3]);
  int n = (int)(w[4]);

  int input_buffer_samples = x->x_input_buffer_samples;
  int input_buffer_mask = input_buffer_samples - 1;
  int write_phase = x->x_write_phase & input_buffer_mask;
  int read_phase = x->x_read_phase;
  t_looper_state state = x->x_state;

  t_sample *vp = x->x_input_buffer;

  while (n--) {
    t_sample f = *in1++;
    t_sample ls = 0.0f;
    int distance_to_end;

    switch (state) {
      case STATE_IDLE:
        vp[write_phase] = f;
        break;
      case STATE_RECORDING:
        vp[write_phase] = f;
        ls = f;
        break;
      case STATE_PLAYING:
        x->x_loop_pos++;
        if (x->x_loop_pos >= x->x_loop_length) x->x_loop_pos = 0;
        ls = vp[read_phase] * x->x_window_buffer[x->x_loop_pos];
        break;
      }

    *out++ = ls;

    write_phase = (write_phase + 1) & input_buffer_mask;
    if (state == STATE_PLAYING) {
      read_phase = (read_phase + 1) & input_buffer_mask;
      int loop_pos = (read_phase - x->x_loop_start) & input_buffer_mask;
      if (loop_pos >= x->x_loop_length) {
        read_phase = x->x_loop_start;
      }
    }
  }

  x->x_read_phase = read_phase;
  x->x_write_phase = write_phase;
  return (w+5);
}

static void looper_dsp(t_looper *x, t_signal **sp)
{
  dsp_add(looper_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_length);
  set_system_params(x, sp[0]->s_length, sp[0]->s_sr);
  input_buffer_update(x);
}

static void looper_bang(t_looper *x)
{
  x->x_state = (x->x_state == STATE_RECORDING) ? STATE_PLAYING : STATE_RECORDING;

  if (x->x_state == STATE_RECORDING) {
    x->x_loop_start = x->x_write_phase;
    x->x_loop_length = 0;
  } else {
    x->x_read_phase = x->x_loop_start;
    x->x_loop_pos = 0;
    x->x_loop_length = (x->x_write_phase - x->x_loop_start) & x->x_input_buffer_samples - 1;
    for (int i = 0; i < x->x_loop_length; i++) {
      float phase = (float)i / (x->x_loop_length - 1);
      x->x_window_buffer[i] = 0.5f * (1.0f - cos(2.0f * M_PI * phase));
    }
  }
}

 static void looper_idle(t_looper *x)
{
  x->x_state = STATE_IDLE;
}


static void looper_free(t_looper *x) {
  if (x->x_input_buffer != NULL) {
    freebytes(x->x_input_buffer, x->x_input_buffer_samples * sizeof(t_sample));
    x->x_input_buffer = NULL;
  }

  if (x->x_window_buffer != NULL) {
    freebytes(x->x_window_buffer, x->x_input_buffer_samples * sizeof(t_sample));
    x->x_window_buffer = NULL;
  }
}

void looper_tilde_setup(void)
{
  looper_class = class_new(gensym("looper~"),
                           (t_newmethod)looper_new,
                           (t_method)looper_free,
                           sizeof(t_looper),
                           CLASS_DEFAULT,
                           A_DEFFLOAT, 0);

  class_addmethod(looper_class, (t_method)looper_dsp,
                  gensym("dsp"), A_CANT, 0);
  class_addmethod(looper_class, (t_method)looper_idle,
                  gensym("idle"), 0);
  class_addbang(looper_class, looper_bang);

  CLASS_MAINSIGNALIN(looper_class, t_looper, x_f);
}

