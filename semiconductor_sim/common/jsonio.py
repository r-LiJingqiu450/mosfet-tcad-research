"""Strict JSON for files and console: one encoder, numpy support, no NaN/Infinity."""
import json
import os
from pathlib import Path
import numpy as np


def json_default(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,Path):return str(value)
    raise TypeError(f'Unsupported JSON type: {type(value).__name__}')


def dumps(value, **kwargs):
    return json.dumps(value,default=json_default,allow_nan=False,ensure_ascii=False,**kwargs)


def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    text=dumps(value,indent=2)+'\n'
    temp=path.with_name(path.name+'.tmp')
    with temp.open('w',encoding='utf-8') as f:
        f.write(text);f.flush();os.fsync(f.fileno())
    temp.replace(path)
