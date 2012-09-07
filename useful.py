def toNum(s):
    try:
	return int(s)
    except ValueError:
	return float(s)


class dotdict(dict):
    def __getattr__(self, name):
	return self[name]
    def __setattr__(self, name, value):
	self[name] = value
	return self
