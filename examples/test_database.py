"""
Open and view items in a tetra3 database.
"""

from numpy import load

database = 'tetra3/default_database.npz'
data = load(database)

lst = data.files 

for item in lst: 
    print(item) 
    print(data[item])
    try:
        print(len(data[item]))
    except:
        pass
