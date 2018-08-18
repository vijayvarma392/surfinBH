import fit_evaluators

FIT_CLASSES = {
    "3dq8": fit_evaluators.Fit3dq8,
    "7dq2": fit_evaluators.Fit7dq2,
    }

def LoadFits(name):
    if name not in FIT_CLASSES.keys():
        raise Exception('Invalid fit name : %s'%name)
    else:
        return FIT_CLASSES[name](name)
