import sqlalchemy as sql
from sqlalchemy.orm import sessionmaker

import models

ENGINE = None
SESSION = None

def connect(filename=None, protocol='sqlite:///'):
    global ENGINE, SESSION
    if filename:
        ENGINE = ENGINE or sql.create_engine(protocol+filename)
        SESSION = SESSION or sessionmaker(bind=ENGINE)
    return SESSION()
