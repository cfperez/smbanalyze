import sqlalchemy as sql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, ForeignKey, Float, LargeBinary
from sqlalchemy.orm import relationship, backref


Base = declarative_base()

class Refold(Base):
  __tablename__ = 'refolds'

  id = Column(Integer, primary_key=True)
  molecule_id = Column(Integer, ForeignKey('molecules.id'))
  refold_time = Column(Float)
  pull = Column(Integer)

  @property
  def construct(self):
    return self.molecule.construct

  trap_data = Column(LargeBinary, nullable=False)
  fret_data = Column(LargeBinary)

#molecule = relationship("Molecule", backref=backref('refolds', order_by=(refold_time,pull)))

  def __repr__(self):
    return "<Refold(construct='{}', refold_time={}, pull={})>".format(
        self.molecule.construct, self.refold_time, self.pull)


class Molecule(Base):
  __tablename__ = 'molecules'

  id = Column(Integer, primary_key=True)
  construct = Column(String)
  slide_id = Column(Integer)
  mol_id = Column(Integer)
  date = Column(String)

  refolds = relationship("Refold", order_by='(Refold.refold_time,Refold.pull)', backref="molecule")

  def __repr__(self):
    return "Molecule(construct='{}', slide_id={}, mol_id={}, date='{}'".format(
          self.construct, self.slide_id, self.mol_id, self.date)

