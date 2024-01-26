#!/usr/bin/env python3
if __name__ == "__main__":
    import sys
    import my_modules
    
    my_modules.helloworld()
    
    arguments = sys.argv[1:]
    
    print(my_modules.add(arguments[0], arguments[1]))
    