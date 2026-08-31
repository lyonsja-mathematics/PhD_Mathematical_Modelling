def non_negative_args(fn):
    def check_non_negative(*args, **kwargs):
        for a in args:
            if isinstance(a, (int, float)):
                if a < 0:
                    raise ValueError("No, that doesn't work. Input argument must be non-negative")
        for k in kwargs.values():
            if isinstance(k, (int, float)):
                if k < 0:
                    raise ValueError("No, that doesn't work. Input keyword argument must be non-negative")
        return fn(*args, **kwargs)
    return check_non_negative
