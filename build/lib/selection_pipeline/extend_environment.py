#
# Extend environment extend the environment variables
# so that the program runs on any linux machine.
#
#
#
import os

def set_environment(environment_variables):
    for environ, value in environment_variables.items():
       os.environ[environ].append(value)

    print("did it works")




