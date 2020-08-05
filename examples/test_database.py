from numpy import load

data = load('tetra3/default_database.npz')
# data = load('stereo/data/tetra3/mag6.7fov7.npz')

lst = data.files 

for item in lst: 
    print(item) 
    print(data[item])
    try:
        print(len(data[item]))
    except:
        pass
