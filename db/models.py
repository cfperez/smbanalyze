import sqlalchemy as sql
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy import Column, Integer, String, ForeignKey, Float, PickleType
from sqlalchemy.orm import relationship, backref, sessionmaker

__all__ = ["Pull", "Refold", "Molecule"]

def to_datetime_str(date_list):
    Y = date_list[0]
    if Y < 2000:
        Y += 2000
    
    return '-'.join(map(lambda s:'{:02d}'.format(s),[Y]+date_list[1:]))

ENGINE = None
SESSION = None

def connect_sql(filename=None, protocol='sqlite:///'):
    global ENGINE, SESSION
    if filename is None and not ENGINE:
        raise ValueError("Must specify database file to connect to")
    if filename:
        ENGINE = ENGINE or sql.create_engine(protocol+filename)
        SESSION = SESSION or sessionmaker(bind=ENGINE)
    return SESSION and SESSION()
Base = declarative_base()

class Molecule(Base):
  __tablename__ = 'molecules'

  id = Column(Integer, primary_key=True)
  construct = Column(String)
  slide_id = Column(Integer)
  mol_id = Column(Integer)
  date = Column(String)

  def __repr__(self):
    return "Molecule(construct='{}', slide_id={}, mol_id={}, date='{}'".format(
          self.construct, self.slide_id, self.mol_id, self.date)


class DataMixIn(object):
  @declared_attr
  def __tablename__(cls):
    return cls.__name__.lower()
  id = Column(Integer, primary_key=True)
  @declared_attr
  def molecule_id(cls):
    return Column(Integer, ForeignKey('molecules.id'))
  pull = Column(Integer)
  trap_data = Column(PickleType, nullable=False)
  fret_data = Column(PickleType)
  meta_data = Column(PickleType)

class Pull(DataMixIn, Base):
  molecule = relationship("Molecule", backref=backref('pulls'))

  def __repr__(self):
    return "<Pull(construct='{}',pull={}>".format(
      self.construct,
      self.pull)

class Refold(DataMixIn, Base):
  refold_time = Column(Float)
  molecule = relationship("Molecule", 
    backref=backref('refolds', order_by='(Refold.refold_time,Refold.pull)'))

  def __repr__(self):
    return "<Refold(construct='{}', refold_time={}, pull={})>".format(
        self.molecule.construct, self.refold_time, self.pull)


def save_refold(exps, conn=session):
    '''Save experiments to database'''
    mol_db = db.Molecule(construct=construct, mol_id=mol_id, 
        slide_id=slide_id, date=date_)
    for p in exps:
        try:
            refold_time = p.metadata['trap.refolding_time']
            pull = p.metadata['trap.current_pull']
            conn.add(db.Refold(molecule=mol_db, refold_time=refold_time, pull=pull,
                     trap_data=p.trap, fret_data=p.fret, meta_data=p.metadata))
        except KeyError:
            print "Expected refolding exp, but got:", p
    conn.commit()
