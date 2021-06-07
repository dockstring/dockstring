from typing import Dict

import pandas as pd

from dockstring import load_target, list_all_target_names
from dockstring.utils import parse_search_box_conf


def compute_volume(config: Dict[str, float]) -> float:
    return config['size_x'] * config['size_y'] * config['size_z']


targets = [load_target(target_name) for target_name in list_all_target_names()]

data = [(target.name, compute_volume(parse_search_box_conf(target.conf_path))) for target in targets]

df = pd.DataFrame(data=data, columns=('target', 'volume'))

df = df.sort_values(by='volume')
print(df)
print()
print(df['volume'].describe())
