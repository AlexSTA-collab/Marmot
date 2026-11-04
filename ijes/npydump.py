"""Output manager 'npydump': the monitor of EdelweissFE, writing the last result of its field output to <label>.npy in
double precision at the end of every increment instead of printing it, and the coordinates of the field's nodes to
<label>_xyz.npy (node fields only). EnSight output is single precision; this is for quantities below its resolution."""
import numpy as np

from edelweissfe.outputmanagers.monitor import OutputManager as Monitor


class OutputManager(Monitor):
    identification = "NpyDump"

    def finalizeIncrement(self, **kwargs):
        for job in self.monitorJobs:
            fo = job["fieldOutput"]
            np.save(f"{self.name}.npy", np.asarray(job["f(x)"](fo.getLastResult()), dtype=np.float64))
            try:
                np.save(f"{self.name}_xyz.npy", np.array([n.coordinates for n in fo.associatedSet], dtype=np.float64))
            except (AttributeError, TypeError):
                pass
