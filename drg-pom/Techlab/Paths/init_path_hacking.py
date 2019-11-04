" init.py path hacking "

import os
import os.path

user_dir = os.path.expanduser('~\.dirforpath')
if os.path.exists(user_dir):
    plugins = os.listdir(user_dir)
    for plugin in plugins:
        __path__.append(os.path.join(user_dir, plugin))