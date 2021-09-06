from .rate import rate
import numpy as np

class DetConfig:
    def __init__(self, B, S=([0.,1.],[1.,1.]), time_window="auto", name=None):
        """
        Detector configuration, containing signal `S` and background `B` event rates.

        Args:
            B(rate): background rate vs. time
            S(rate): expected signal event rate vs. time
            time_window ("auto" or tuple(float, float)):
                The limits around t0 in which to take the signal.
                If "auto", then try to take the full range from `S` (via S.range)
            name(str): optional name for the detector configuration
        """
        self.S = rate(S)
        if(time_window=="auto"):
            time_window = self.S.range

        if np.any(np.isinf(time_window)):
            raise ValueError(f'Cannot work with infinite time window: {time_window}')

        self.S0=self.S.integral(*time_window)
        self.B = rate(B) 

        self.time_window=np.array(time_window)
        self.name=name
    
    def __str__(self):
        tw = self.time_window
        Ns = self.S.integral(*tw)
        Nb = self.B.integral(*tw)
        return f"<DetConfig '{self.name}' = "+f"time_window=({tw[0]}, {tw[1]}), "+\
                f"S/B = {Ns}/{Nb} >"
