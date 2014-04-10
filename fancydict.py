class dotdict(dict):
  def __init__(self,*args,**kwargs):
    super(dotdict,self).__init__(*args,**kwargs)

  def __setitem__(self, key, val):
    if not str(key).replace('_','').isalnum():
      raise KeyError, "Key must be alphanumeric"
    super(dotdict,self).__setitem__(key,val)

  def __getitem__(self,key):
    key,_,the_rest = key.partition('.')
    if not self.has_key(key) and \
       not str(key).startswith('_'):
      self[key] = dotdict()
    if the_rest:
      return super(dotdict,self).__getitem__(key)[the_rest]
    return super(dotdict,self).__getitem__(key)

  def __getattr__(self, name):
    if not self.__dict__.has_key(name):
      try:
        return self.__getitem__(name)
      except KeyError:
        raise AttributeError('Dotdict has no attribute %s' % name)
    else: return super(dotdict,self).__getattr__(name)

  def __setattr__(self, name, value):
    if name.startswith('_'):
      super(dotdict,self).__setattr__(name,value)
    else:
      self[name] = value

class nesteddict(dict):
  @classmethod
  def from_dict(cls, *args, **kwargs):
    dict_ = dict(*args, **kwargs)
    new_dict = cls()
    for key,val in dict_.iteritems():
      k,_,more = key.partition('.')
      context = new_dict.setdefault(k, cls()) if more else new_dict
      if more:
        k = more
      context[k] = cls.from_dict(val) if isinstance(val, dict) else val
    return new_dict

  def __getitem__(self, key):
    key,_,more = key.partition('.')
    val = super(nesteddict,self).__getitem__(key)
    if more:
      if isinstance(val,dict):
        return val[more]
      raise KeyError('"{}" not found in key {}'.format(more,key))
    return val

  def __setitem__(self, key, val):
    key,_,more = key.partition('.')
    if more:
      context = self[key]
      if isinstance(context,dict):
        context[more] = val
        return
    super(nesteddict,self).__setitem__(key, val) 