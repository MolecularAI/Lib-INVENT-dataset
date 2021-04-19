#!/usr/bin/env python

import sys
import json

from dacite import from_dict

from dto import GeneralConfiguration
from manager import Manager

if __name__ == "__main__":

    with open(sys.argv[1]) as f:
        json_input = f.read().replace('\r', '').replace('\n', '')

    configuration = {}
    try:
        configuration = json.loads(json_input)
    except (ValueError, KeyError, TypeError):
        print("JSON format error")
    else:
        general_config = from_dict(data_class=GeneralConfiguration, data=configuration)
        manager = Manager(general_config)
        manager.run()
