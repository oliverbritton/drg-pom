COMMENT
iramp.mod

Delivers a ramp current that starts at t = delay >=0.
The current starts at zero and increases linearly until t = delay+dur.
The current maximum amplitude is given by amp.

Uses event delivery system to ensure compatibility with adaptive integration.

ENDCOMMENT

NEURON {
  POINT_PROCESS IRamp
  RANGE delay, dur, amp
  ELECTRODE_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  delay (ms)
  dur (ms)
  amp (nA)
}

ASSIGNED {
  i (nA)
  on (1)
}

INITIAL {
  i = 0
  on = 0

  if (delay<0) { delay=0 }
  if (dur<0) { dur=0 }

  : do nothing if dur == 0
  if (dur>0) {
    net_send(delay, 1)  : to turn it on and start ramp
  }
}


BREAKPOINT {
  if (on==0) {
    i = 0
  } else {
    i = amp * (t - delay)/dur
  }
}

NET_RECEIVE (w) {
  : respond only to self-events with flag = 1
  if (flag == 1) {
    if (on==0) {
      on = 1  : turn it on
      net_send(dur, 1)  : to stop frequency ramp, freezing frequency at f1
    } else {
      on = 0  : turn it off
    }
  }
}
