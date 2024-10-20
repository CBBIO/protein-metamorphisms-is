from datetime import datetime

from pgvector.sqlalchemy import Vector
from sqlalchemy import (Column, Integer, String, Date, ForeignKey, DateTime,
                        func, Float, Boolean, Index, UniqueConstraint)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import relationship, declarative_base, mapped_column

Base = declarative_base()