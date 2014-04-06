'''datetime functions/wrappers

mongoDB only handles (natively) datetime objects, not date object. This
is an issue when trying to find things by only YMD (datetime objects
include HMSuS.)
'''

import datetime

def today():
  '''Return datetime object with YMD and HMSuS=0''' 
  d = datetime.date.today()
  return date(d.year,d.month,d.day)

def date(y,m,d,h=0,mn=0,s=0):
  '''Return datetime object with default HMS=0
  Year "y" as two digits (e.g. 14) is interpreted as 2014'''
  if y < 2000:
      y += 2000
  return datetime.datetime(y,m,d,h,mn,s)

def to_date(date_time):
  return date(date_time.year, date_time.month, date_time.day)
  
def str_to_date(date_string):
  date_num = map(int, date_string.split('.'))
  if date_num[0] < 2000:
      date_num[0] += 2000
  return datetime.date(*date_num)