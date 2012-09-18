def toNum(s):
  if s is None:
	return None
  try:
	return int(s)
  except ValueError:
	return float(s)

def toInt(s):
	s = toNum(s)
	return int(s) if s else None

class dotdict(dict):
  def __setitem__(self, key, val):
	if not str(key).replace('_','').isalnum():
	  raise KeyError, "Key must be alphanumeric"
	super(dotdict,self).__setitem__(key,val)

  def __getitem__(self,key):
	if not self.has_key(key):
	  self[key] = dotdict()
	return super(dotdict,self).__getitem__(key)

  def __getattr__(self, name):
	return self[name]

  def __setattr__(self, name, value):
	self[name] = value
	return self
