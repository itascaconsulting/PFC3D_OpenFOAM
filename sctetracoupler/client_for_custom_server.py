import customsocket
import create_mesh
#reload(customsocket)


#client = customsocket.CustomSocketClient(tcp_id='192.168.1.11',tcp_port=3333)
client = customsocket.CustomSocketClient(tcp_id='127.0.0.1',tcp_port=3333)

client.start()


# UPDATE MESH
# receive number of nodes
nbNodes = client.read_data()
print "Number of nodes: ",nbNodes

# receive node information
nodes = {}
for n in range(0,nbNodes):
  i = client.read_data() # node id
  x = client.read_data() # node x 
  y = client.read_data() # node y
  z = client.read_data() # node z
  nodes[i] = [x,y,z]
  #print "coord node {} : {} {} {}".format(n,x,y,z)

# receive connectivity information
nbElem = client.read_data()
print "Number of elements: ",nbElem
elements_nodes = {}
for n in range(0,nbElem):
  i  = client.read_data()  # element ID (should be n)
  nn = client.read_data()  # total number of nodes
  idx = []
  for j in range(0,nn):
    nj = client.read_data()  # index of node j
    idx.append(nj)
  elements_nodes[i] = idx
#  print "connection {} : {}".format(n,ic)


# UPDATE MESHELEM
# receive number of elements
nbElem = client.read_data()
print "Number of elements: ",nbElem
elements_geom = {}
for n in range(0,nbElem):
  x = client.read_data()
  y = client.read_data()
  z = client.read_data()
  v = client.read_data()
  elements_geom[n] = [x,y,z,v]
  #print "element {} : x = {} y = {} z = {} - vol = {}".format(n,x,y,z,v)


# UPDATE FIELD
# receive number of elements
nbElem = client.read_data()
print "Number of elements: ",nbElem
elements_field = {}
for n in range(0,nbElem):
  u = client.read_data()
  v = client.read_data()
  w = client.read_data()
  p = client.read_data()
  d = client.read_data()
  vi = client.read_data()
  elements_field[n] = [u,v,w,p,d,vi]
  #print "element {} : u = {} v = {} w = {} p = {} d = {} vi = {}".format(n,u,v,w,p,d,vi)

#send communication status (success)
client.send_data(1)

client.close()

# create mesh as a collection of geometry sets for visualization
create_mesh.create_mesh(nodes,elements_nodes)