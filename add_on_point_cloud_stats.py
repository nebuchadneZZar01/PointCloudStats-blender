# author: Michele Ferro
# mail: michele.ferro1998@libero.it
# github: github.com/nebuchadneZZar01

bl_info = {
    "name": "Point Cloud Stats",
    "author": "nebuchadnezzar",
    "version": (1, 5),
    "blender": (2, 90, 0),
    "location": "Scene > Point Cloud Stats",
    "description": "Retrieves meshes informations, saves them in a CSV file and creates a collection of 'Pills' that summarizes them.",
    "warning": "",
    "doc_url": "",
    "category": "Add Mesh",
}

# Importazioni delle librerie necessarie
# mathutils contiene utili metodi di conversione delle matrici (da rotazione ad angolo di eulero)
# operator contiene utili operatori add, abs, ecc.
# itertools contiene utili operatori di ritaglio dei dizionari

import bpy
import math
import mathutils 
import itertools
import operator
import csv

#======================================================================
#                  FUNZIONI DI CALCOLO
#======================================================================

# --- FUNZIONI DI UTILITÀ GENERALE ---
def distEuclide(vert_source, vert_dest):
    xS = vert_source[0]
    yS = vert_source[1]
    zS = vert_source[2]

    xD = vert_dest[0]
    yD = vert_dest[1]
    zD = vert_dest[2]

    return abs(((xD-xS)**2 + (yD-yS)**2 + (zD-zS)**2)**1/2)

def isNegative(num):
    if num < 0: return True
    else: return False

def isPositive(num):
    if num >= 0: return True
    else: return False

# --- INDICI DI TENDENZA CENTRALE ---

# calcolo media elementi di una lista
def media(L):
    if len(L) != 0: return sum(L)/len(L)
    else: return 0.0

# calcolo varianza elementi di una lista
def stDev(L):
    m = media(L)
    L = [(x-m)*(x-m) for x in L]
    return pow(sum(L)/len(L),1/2)


# --- FUNZIONI DI CALCOLO DEGLI AUTOVALORI/VETTORI FORNITE DAL PROF. GALLO ---
# A sia una matrice 3 x 3
# descritta da una lista di 3 liste
# ciascuna delle quali si compone di 3 elementi

def modulo(v):
    w = [x*x for x in v]
    return (sum(w))**(1/2)

def prodotto(A,B):
    # A e B sono due liste di liste
    # che descrivono due matrici (MxN) e (NxK)
    rowA = len(A)                 # righe di A
    colA = len(A[0])              # colonne di A
    rowB = len(B)                 # righe di B
    colB = len(B[0])              # colonne di B
    C = []
    for i in range(rowA):
        C.append([0 for x in range(colB)])
         # calcoliamo i prodotti 
    for i in range(rowA):
        for j in range(colB):
            for k in range(rowA):
                C[i][j] = C[i][j] + A[i][k]*B[k][j]
    # restituiamo C
    return C
     
def decomposizioneQR(A):
    # deve restituire due matrici Q e R
    # Q ortonormale e R upper triangolare
    # CALCOLO prima colonna di Q
    # è eguale alla prima colonna di A ma normalizzata
    m1 = modulo([A[0][0], A[1][0], A[2][0]])
    q00 = A[0][0]/m1
    q10 = A[1][0]/m1
    q20 = A[2][0]/m1
    # CALCOLO seconda colonna di Q
    # inizio con il calcolare quanto è lunga la proiezione della
    # seconda colonna di A sulla prima colonna di Q
    s = A[0][1]*q00 + A[1][1]*q10 + A[2][1]*q20
    # sottraggo dalla seconda colonna di A un vettore parallelo a 
    # alla prima colonna di Q e lungo s, quel che
    # rimane è perpendicolare alla prima colonna di Q,e lo normalizzo
    a01 = A[0][1] - s*q00
    a11 = A[1][1]-s*q10
    a21 = A[2][1]-s*q20
    m2 = modulo([a01,a11,a21])
    q01 = a01/m2
    q11 = a11/m2
    q21 = a21/m2
    # CALCOLO terza colonna di Q
    # calcolo la lunghezza della proiezione della terza colonna di A sulla 
    # prima di Q
    s0 = A[0][2]*q00 + A[1][2]*q10 + A[2][2]*q20
    # calcolo la lunghezza della terza colonna di A sulla seconda di Q
    s1 = A[0][2]*q01 + A[1][2]*q11 + A[2][2]*q21
    # sottraggo le componenti così trovate e trovo il terzo vettore che normalizzo
    a02 = A[0][2] - s0*q00 - s1*q01
    a12 = A[1][2] - s0*q10 - s1*q11
    a22 = A[2][2] - s0*q20 - s1*q21
    m3 = modulo([a02, a12, a22])
    q02 = a02/m3
    q12 = a12/m3
    q22 = a22/m3
    # calcolo matrice R
    r10 = 0
    r20 = 0
    r21 = 0
    r00 = m1
    r11 = m2
    r22 = m3
    r01 = s
    r02 = s0
    r12 = s1
    # costruzione e finale ?????????  occhio agli indici
    Q=[[q00,q01,q02], [q10,q11,q12], [q20,q21,q22]]
    R=[[r00,r01,r02], [r10,r11,r12], [r20,r21,r22]]
    return [Q,R]

def utu(M):
    # calcola approximate eingenvalue e eigenvector
    # solitamente k = 10, proviamo con valori più bassi
    U = [[1,0,0], [0,1,0], [0,0,1]]
    A = M.copy()
    k = 5
    for i in range(k):
        qr = decomposizioneQR(A)
        A = prodotto(qr[1],qr[0])
        U1 = U.copy()
        U = prodotto(U1,qr[0])
    return (A,U)

# Usata nella normalizzazione degli autovalori
def normalize(eig):
    tot = sum([abs(x) for x in eig])
    eig = [x/tot for x in eig]
    return eig


# --- FUNZIONI DI CALCOLO DEGLI ELEMENTI DEL TENSORE DI INERZIA --- 

def momentoInerziaDiagonale(V,W):
    V = [x*x for x in V]
    W = [x*x for x in W]
    return sum(V)+sum(W)
    
def momentoInerziaMisto(V,W):
    inerzia = 0
    for i in range(len(V)):
        inerzia = inerzia - (V[i]*W[i])
    return inerzia



#======================================================================
#                  FUNZIONI DI COLORAZIONE
#======================================================================

# --- GENERAZIONE DEI MATERIALI ---
# Crea preventivamente cinque materiali, ognuno con un suo valore RGB
# e con roughness massima, per poi inserirli in una lista
def generateMaterial():
    mat1 = bpy.data.materials.new(name = "000_020_MAT")
    mat2 = bpy.data.materials.new(name = "020_040_MAT")
    mat3 = bpy.data.materials.new(name = "040_060_MAT")
    mat4 = bpy.data.materials.new(name = "060_080_MAT")
    mat5 = bpy.data.materials.new(name = "080_100_MAT")

    mat1.roughness = mat2.roughness = mat3.roughness = mat4.roughness = mat5.roughness = 1

    mat1.diffuse_color = (0.8,0.8,0.9,1)
    mat2.diffuse_color = (0.6,0.7,0.8,1)
    mat3.diffuse_color = (0.3,0.5,0.7,1)
    mat4.diffuse_color = (0.0,0.3,0.4,1)
    mat5.diffuse_color = (0.0,0.1,0.5,1)

    matList = [mat1, mat2, mat3, mat4, mat5]

    return matList 

# --- ASSEGNAZIONE DEI MATERIALI ---
# preso in input una pillola, un parametro e una lista di materiali
# viene assegnato un materiale in base al range di appartenenza di tale parametro in modo lineare.
# Al momento tale funzione è configurata per funzionare sulla base del prodotto degli scarti quadratici
# delle normali, il cui massimo valore è risultato essere circa 0.20, pertanto è stato suddiviso
# in modo lineare l'intervallo [0.0,0.20] in cinque sottoclassi.
def assignMaterialLinear(pillObj, parameter, materials):   
    if parameter >= 0.00 and parameter < 0.20:
        pillObj.data.materials.append(materials[0])
    elif parameter >= 0.20 and parameter < 0.40:
        pillObj.data.materials.append(materials[1])
    elif parameter >= 0.40 and parameter < 0.60:
        pillObj.data.materials.append(materials[2])
    elif parameter >= 0.60 and parameter < 0.80:
        pillObj.data.materials.append(materials[3])
    elif parameter >= 0.80 and parameter <= 1:
        pillObj.data.materials.append(materials[4])


#======================================================================
# GESTIONE DEL FILE CSV
# Funzione che gestisce la creazione del file CSV
csvFile = None
writer = None
meshesData = dict()

def createCSV(filepath): 
    global csvFile
    global writer
    
    csvFile = open(filepath, 'w', newline='')

    # compatibilità con tabulazione di excel per dati separati da virgola
    csvFile.write("sep=,\n")

    # lista dei campi del foglio csv
    fieldnames = ['object', 'collection', 'nVert', 'nFace',\
                  'totSurface', 'Vol', 'planarity',
                  'mediaVX', 'mediaVY', 'mediaVZ',\
                  'varVX', 'varVY', 'varVZ',\
                  'mediaNX', 'mediaNY', 'mediaNZ',\
                  'varNX', 'varNY', 'varNZ',\
                  'eigv1', 'eigv2', 'eigv3',\
                  'x1', 'y1', 'z1',\
                  'x2', 'y2', 'z2',\
                  'x3', 'y3', 'z3',\
                  'XN_YN_ZN_x', 'XN_YN_ZN_y', 'XN_YN_ZN_z',\
                  'XN_YN_ZP_x', 'XN_YN_ZP_y', 'XN_YN_ZP_z',\
                  'XN_YP_ZN_x', 'XN_YP_ZN_y', 'XN_YP_ZN_z',\
                  'XN_YP_ZP_x', 'XN_YP_ZP_y', 'XN_YP_ZP_z',\
                  'XP_YN_ZN_x', 'XP_YN_ZN_y', 'XP_YN_ZN_z',\
                  'XP_YN_ZP_x', 'XP_YN_ZP_y', 'XP_YN_ZP_z',\
                  'XP_YP_ZN_x', 'XP_YP_ZN_y', 'XP_YP_ZN_z',\
                  'XP_YP_ZP_x', 'XP_YP_ZP_y', 'XP_YP_ZP_z',\
                  'objType'] 
    
    writer = csv.DictWriter(csvFile, fieldnames=fieldnames)
    writer.writeheader()

    return csvFile, writer

def writeCSV(filepath):
    csvFile, writer = createCSV(filepath)

    try:
        for obj in meshesData:
            writer.writerow(meshesData[obj])
    except:
        print("Error in writer.writerow(): csv file was not created!")

    try:
        csvFile.close()
    except:
        print("Error in csv.close(): csv file was not created!")
    
#========================================================

#======================================================================
#               SCRIPT PRINCIPALE da eseguire
#        (tutte le funzioni create sopra servono qui!)
#======================================================================
def main(context):
    global csvFile
    global writer
    global meshesData

    #========================================================
    # GESTIONE CREAZIONE DELLA COLLEZIONE DI PILLOLE
    # Nel caso in cui l'utente volesse ri-eseguire il calcolo, 
    # la vecchia collezione verrà cancellata e sostituita con
    # una nuova collezione di pillole.
    cols = bpy.data.collections
    for col in cols:
        if col.name == "PillsCollection":
            bpy.data.collections.remove(col)

    myCol = bpy.data.collections.new("PillsCollection")
    bpy.context.scene.collection.children.link(myCol)
    #========================================================
       

    #========================================================
    # CALCOLI!
    # carico nella lista objs tutte le mesh visibili presenti nella scena     
    objs = bpy.context.visible_objects.copy()

    # creo i materiali che serviranno ad analizzare la rugosità delle pietre
    pillsColors = generateMaterial()

    # I calcoli vengono eseguiti per ogni oggetto visibile in scena
    # Inizialmente è necessario assicurarsi che i cambiamenti effettuati in
    # scena siano applicati, così che il data object della mesh sia aggiornato
    for obj in objs:
        if obj.type == "MESH":
            bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
            colName = obj.users_collection[0].name_full
            obj['Object Type'] = ""
            obj_data = obj.data
            
            # --- SEZIONE RELATIVA AI VERTICI ---
            # prendo il riferimento dei vertici del data object in viewport ed il loro numero
            """Sono necesarie non le coordinate locali, bensì quelle globali. 
            Per trovarle devo moltiplicare ogni vettore che rappresenta il vertice 
            con la matrice di trasformazione del mondo di blender
            translated_vertex = obj.matrix_world @ obj_verts[i].co
            Questa operazione deve essere fatto per ogni vertice presente nella mesh, 
            pertento è necessario trasferirte tutti i vertici all'interno 
            di una nuova struttura dati obj_verts e poi effettuare
            quanto già fatto in precedenza nello script"""
            
            obj_verts = obj_data.vertices
            m_world = obj.matrix_world              # matrice di trasformazione globale
            #m_i_world = m_world.inverted()          # inversa della matrice di trasformazione globale
            nVerts = len(obj_verts)

            x_verts = []
            y_verts = []
            z_verts = []

            x_vertsAvg = y_vertsAvg = z_vertsAvg = 0
            x_vertsStDev = y_vertsStDev = z_vertsStDev = 0

            # raccolgo le coordinate GLOBALI: il passaggio di sistema di riferimento viene eseguito moltiplicando i vettori che rappresentano 
            # i vertici per la matrice di trasformazione del mondo di blender dei vertici e le inserisco in tre liste, una per ogni asse
            for vert in obj_verts:
                w_vert = m_world @ vert.co
                x_verts.append(w_vert.x)
                y_verts.append(w_vert.y)
                z_verts.append(w_vert.z)

            # CALCOLO DELLE MEDIE PER LE COORDINATE IN X, Y, Z (Il punto corrisponderà poi al baricentro GEOMETRICO della mesh)
            x_vertsAvg = media(x_verts)
            y_vertsAvg = media(y_verts)
            z_vertsAvg = media(z_verts)

            vertsAvg = (x_vertsAvg, y_vertsAvg, z_vertsAvg)

            # CALCOLO DELLE DEVIAZIONI STANDARD PER LE COORDINATE IN X, Y, Z
            x_vertsStDev = stDev(x_verts)
            y_vertsStDev = stDev(y_verts)
            z_vertsStDev = stDev(z_verts)

            vertsStDev = (x_vertsStDev, y_vertsStDev, z_vertsStDev)
            volume = x_vertsStDev * y_vertsStDev * z_vertsStDev
            rSphere = volume**(1/3)

            # --- SEZIONE RELATIVA ALLE NORMALI ---
            # prendo il riferimento delle facce del data object (la mesh) in viewport
            obj_polygons = obj_data.polygons
            nFaces = len(obj_polygons)

            x_normals = []
            y_normals = []
            z_normals = []
            x_normalsAvg = y_normalsAvg = z_normalsAvg = 0
            totSurface = 0
            
            # calcolo la superficie totale della mesh
            for p in obj_polygons:
                totSurface += p.area

            # inserisco le normali in delle liste, una per ogni coordinata
            for p in obj_polygons:
                normals = p.normal
                x_normals.append(normals[0])
                y_normals.append(normals[1])
                z_normals.append(normals[2])

            # CALCOLO DELLE MEDIE PER LE COORDINATE x, y e z delle normali
            x_normalsAvg = media(x_normals)
            y_normalsAvg = media(y_normals)
            z_normalsAvg = media(z_normals)

            normalsAvg = (x_normalsAvg, y_normalsAvg, z_normalsAvg)

            # CALCOLO DELLE DEVIAZIONI STANDARD x, y e z DELLE NORMALI
            x_normalsStDev = stDev(x_normals)
            y_normalsStDev = stDev(y_normals)
            z_normalsStDev = stDev(z_normals)

            normals_StDev = (x_normalsStDev, y_normalsStDev, z_normalsStDev)
            rough = x_normalsStDev * y_normalsStDev * z_normalsStDev
            
            # --- CREAZIONE ISTOGRAMMA NORMALI LOCALI ---
            """L'istogramma di una pietra consiste in uno schema di 24 valori totali, nello specifico otto bin con tre valori x, y, z. Ad ogni bin, corrisponde un quadrante nello spazio. Vengono quindi analizzate
            le componenti alle normali delle superfici di ogni mesh: per ogni normale, viene esaminato il quadrante di appartenenza, e le singole componenti vengono inserite in delle liste (una per coordinata
            appartenente al rispettivo quadrante). Per indicare i nomi dei vari quadranti, è stato scelto un approccio secondo cui vengono indicati gli assi e i versi che descrivono lo specifico quadrante.
            Viene poi calcolata la normale media in ogni quadrante."""

            # vettori contenenti le normali sugli ottanti per ogni coordinata
            XN_YN_ZN_x = []
            XN_YN_ZN_y = []    
            XN_YN_ZN_z = []

            XN_YN_ZP_x = []    
            XN_YN_ZP_y = []    
            XN_YN_ZP_z = []

            XN_YP_ZN_x = []
            XN_YP_ZN_y = []
            XN_YP_ZN_z = []

            XN_YP_ZP_x = []
            XN_YP_ZP_y = []
            XN_YP_ZP_z = []

            XP_YN_ZN_x = []
            XP_YN_ZN_y = []
            XP_YN_ZN_z = []

            XP_YN_ZP_x = []
            XP_YN_ZP_y = []
            XP_YN_ZP_z = []

            XP_YP_ZN_x = []
            XP_YP_ZN_y = []
            XP_YP_ZN_z = []

            XP_YP_ZP_x = []
            XP_YP_ZP_y = []
            XP_YP_ZP_z = []

            for p in obj_polygons:
                if isNegative(p.normal[0]) and isNegative(p.normal[1]) and isNegative(p.normal[2]):
                    XN_YN_ZN_x.append(p.normal[0])
                    XN_YN_ZN_y.append(p.normal[1])
                    XN_YN_ZN_z.append(p.normal[2])
                elif isNegative(p.normal[0]) and isNegative(p.normal[1]) and isPositive(p.normal[2]):
                    XN_YN_ZP_x.append(p.normal[0])
                    XN_YN_ZP_y.append(p.normal[1])
                    XN_YN_ZP_z.append(p.normal[2])
                elif isNegative(p.normal[0]) and isPositive(p.normal[1]) and isNegative(p.normal[2]):
                    XN_YP_ZN_x.append(p.normal[0])
                    XN_YP_ZN_y.append(p.normal[1])
                    XN_YP_ZN_z.append(p.normal[2])
                elif isNegative(p.normal[0]) and isPositive(p.normal[1]) and isPositive(p.normal[2]):
                    XN_YP_ZP_x.append(p.normal[0])
                    XN_YP_ZP_y.append(p.normal[1])
                    XN_YP_ZP_z.append(p.normal[2])
                elif isPositive(p.normal[0]) and isNegative(p.normal[1]) and isNegative(p.normal[2]):
                    XP_YN_ZN_x.append(p.normal[0])
                    XP_YN_ZN_y.append(p.normal[1])
                    XP_YN_ZN_z.append(p.normal[2])
                elif isPositive(p.normal[0]) and isNegative(p.normal[1]) and isPositive(p.normal[2]):
                    XP_YN_ZP_x.append(p.normal[0])
                    XP_YN_ZP_y.append(p.normal[1])
                    XP_YN_ZP_z.append(p.normal[2])
                elif isPositive(p.normal[0]) and isPositive(p.normal[1]) and isNegative(p.normal[2]):
                    XP_YP_ZN_x.append(p.normal[0])
                    XP_YP_ZN_y.append(p.normal[1])
                    XP_YP_ZN_z.append(p.normal[2])
                elif isPositive(p.normal[0]) and isPositive(p.normal[1]) and isPositive(p.normal[2]):
                    XP_YP_ZP_x.append(p.normal[0])
                    XP_YP_ZP_y.append(p.normal[1])
                    XP_YP_ZP_z.append(p.normal[2])

            # normale media su ogni ottante per coordinata
            avg_XN_YN_ZN_x = media(XN_YN_ZN_x)
            avg_XN_YN_ZN_y = media(XN_YN_ZN_y) 
            avg_XN_YN_ZN_z = media(XN_YN_ZN_z)

            avg_XN_YN_ZP_x = media(XN_YN_ZP_x)
            avg_XN_YN_ZP_y = media(XN_YN_ZP_y) 
            avg_XN_YN_ZP_z = media(XN_YN_ZP_z)

            avg_XN_YP_ZN_x = media(XN_YP_ZN_x)
            avg_XN_YP_ZN_y = media(XN_YP_ZN_y) 
            avg_XN_YP_ZN_z = media(XN_YP_ZN_z)

            avg_XN_YP_ZP_x = media(XN_YP_ZP_x)
            avg_XN_YP_ZP_y = media(XN_YP_ZP_y)
            avg_XN_YP_ZP_z = media(XN_YP_ZP_z)

            avg_XP_YN_ZN_x = media(XP_YN_ZN_x)
            avg_XP_YN_ZN_y = media(XP_YN_ZN_y) 
            avg_XP_YN_ZN_z = media(XP_YN_ZN_z)

            avg_XP_YN_ZP_x = media(XP_YN_ZP_x)
            avg_XP_YN_ZP_y = media(XP_YN_ZP_y) 
            avg_XP_YN_ZP_z = media(XP_YN_ZP_z)

            avg_XP_YP_ZN_x = media(XP_YP_ZN_x)
            avg_XP_YP_ZN_y = media(XP_YP_ZN_y) 
            avg_XP_YP_ZN_z = media(XP_YP_ZN_z)

            avg_XP_YP_ZP_x = media(XP_YP_ZP_x)
            avg_XP_YP_ZP_y = media(XP_YP_ZP_y) 
            avg_XP_YP_ZP_z = media(XP_YP_ZP_z)

            # --- CALCOLO DEL TENSORE DI INERZIA ---
            # tolgo l'offset del centro dell'oggetto
            x_v=[x-x_vertsAvg for x in x_verts]
            y_v=[y-y_vertsAvg for y in y_verts]
            z_v=[z-z_vertsAvg for z in z_verts]
               
            # stima raggio della sfera minima centrata nel baricentro che include tutti i vertici della mesh
            raggi = []
            for i in range(len(x_v)):
                raggi.append(x_v[i]**2 + y_v[i]**2 + z_v[i]**2)
            raggioSfera = max(raggi)**(1/2)
            
            I_xx = momentoInerziaDiagonale(y_v, z_v)
            I_yy = momentoInerziaDiagonale(x_v, z_v)
            I_zz = momentoInerziaDiagonale(x_v, y_v)
            I_xy = momentoInerziaMisto(x_v, y_v)
            I_xz = momentoInerziaMisto(x_v, z_v)
            I_yz = momentoInerziaMisto(y_v, z_v)

            T_matrix = [[I_xx, I_xy, I_xz],\
                        [I_xy, I_yy, I_yz],\
                        [I_xz, I_yz, I_zz]]
            
            eigvals, eigvecs = utu(T_matrix)
            autovalori = [eigvals[0][0], eigvals[1][1], eigvals[2][2]]
            
            # --- CALCOLO DELL'INDICE DI PLANARITÀ ---
            """Questo tentativo di implementazione dell'indice di planarità si basa sul seguente approccio: viene presa come riferimento la normale media globale della mesh. Viene generato un piano ad essa
            perpendicolare e successivamente viene calcolata la distanza euclidea tra tutti i punti della mesh e il centro del piano. Tale valore dovrebbe teoricamente corrispondere ad un approssimativo indice
            di planarità della pietra.
            Tuttavia è necessaria una normalizzazione: come indice di normalizzazione è stato scelto il raggio della sfera minima includente la mesh."""

            bpy.ops.mesh.primitive_plane_add(size = 1.0, location = vertsAvg)
            bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')
            plane = bpy.context.active_object
            plane.hide_set(True)
            plane.name = obj.name + '_plane'
            plane.rotation_euler[0] = math.radians(90)
            plane.location = plane.location + mathutils.Vector(normalsAvg)

            dists = []
            for vert in obj_verts:
                w_vert = m_world @ vert.co
                dist = distEuclide(w_vert, plane.location)
                dists.append(dist)

            avgDist = media(dists)
            planarity = avgDist / raggioSfera
            
            # cancello il piano che non dovrò più utilizzare
            bpy.data.objects.remove(plane, do_unlink=True)

            # --- RIEMPIMENTO DIZIONARIO ---
            meshInfo = {'object': obj.name, 'collection': colName, 'nVert': nVerts, 'nFace': nFaces, 
                        'totSurface': round(totSurface,3), 'Vol': round(volume,3), 'planarity': round(planarity,3),\
                        'mediaVX': round(x_vertsAvg,3), 'mediaVY': round(y_vertsAvg,3), 'mediaVZ': round(z_vertsAvg,3),\
                        'varVX': round(x_vertsStDev,3), 'varVY': round(y_vertsStDev,3), 'varVZ': round(z_vertsStDev,3),\
                        'mediaNX': round(x_normalsAvg,3), 'mediaNY': round(y_normalsAvg,3), 'mediaNZ': round(z_normalsAvg,3),\
                        'varNX': round(x_normalsStDev,3), 'varNY': round(y_normalsStDev,3), 'varNZ': round(z_normalsStDev,3),\
                        'eigv1': round(autovalori[0],3), 'eigv2': round(autovalori[1],3), 'eigv3': round(autovalori[2],3),\
                        'x1': round(eigvecs[0][0],3), 'y1': round(eigvecs[0][1],3), 'z1': round(eigvecs[0][2],3),\
                        'x2': round(eigvecs[1][0],3), 'y2': round(eigvecs[1][1],3), 'z2': round(eigvecs[1][2],3),\
                        'x3': round(eigvecs[2][0],3), 'y3': round(eigvecs[2][1],3), 'z3': round(eigvecs[2][2],3),\
                        'XN_YN_ZN_x': abs(round(avg_XN_YN_ZN_x,3)), 'XN_YN_ZN_y': abs(round(avg_XN_YN_ZN_y,3)), 'XN_YN_ZN_z': abs(round(avg_XN_YN_ZN_z,3)),\
                        'XN_YN_ZP_x': abs(round(avg_XN_YN_ZP_x,3)), 'XN_YN_ZP_y': abs(round(avg_XN_YN_ZP_y,3)), 'XN_YN_ZP_z': abs(round(avg_XN_YN_ZP_z,3)),\
                        'XN_YP_ZN_x': abs(round(avg_XN_YP_ZN_x,3)), 'XN_YP_ZN_y': abs(round(avg_XN_YP_ZN_y,3)), 'XN_YP_ZN_z': abs(round(avg_XN_YP_ZN_z,3)),\
                        'XN_YP_ZP_x': abs(round(avg_XN_YP_ZP_x,3)), 'XN_YP_ZP_y': abs(round(avg_XN_YP_ZP_y,3)), 'XN_YP_ZP_z': abs(round(avg_XN_YP_ZP_z,3)),\
                        'XP_YN_ZN_x': abs(round(avg_XP_YN_ZN_x,3)), 'XP_YN_ZN_y': abs(round(avg_XP_YN_ZN_y,3)), 'XP_YN_ZN_z': abs(round(avg_XP_YN_ZN_z,3)),\
                        'XP_YN_ZP_x': abs(round(avg_XP_YN_ZP_x,3)), 'XP_YN_ZP_y': abs(round(avg_XP_YN_ZP_y,3)), 'XP_YN_ZP_z': abs(round(avg_XP_YN_ZP_z,3)),\
                        'XP_YP_ZN_x': abs(round(avg_XP_YP_ZN_x,3)), 'XP_YP_ZN_y': abs(round(avg_XP_YP_ZN_y,3)), 'XP_YP_ZN_z': abs(round(avg_XP_YP_ZN_z,3)),\
                        'XP_YP_ZP_x': abs(round(avg_XP_YP_ZP_x,3)), 'XP_YP_ZP_y': abs(round(avg_XP_YP_ZP_y,3)), 'XP_YP_ZP_z': abs(round(avg_XP_YP_ZP_z,3)),\
                        'objType': obj['Object Type']}

            meshesData[obj.name] = meshInfo

            # --- VISUALIZZAZIONE GRAFICA DATI RACCOLTI ---
            """Uso come sfera da deformare la sfera di raggio minimo centrata nel baricentro
            che contiene la mesh e un fattore di scala "globale" scaleFac
            La normalizzazione va effettuata dopo il calcolo: effettuandola prima si avrebbe una dispersione dei valori, causando l'ingigantimento degli elissoidi.
            Le pillole possono essere colarate: al momento tale feature è impostata usando come parametro il prodotto degli scarti quadratici delle normali, ma è temporaneo, per
            mostrare che è una feature utilizzabile con parametri adeguati."""

            scaleFac = 3
            diags = [(abs(autovalori[0]))**(-1/2), (abs(autovalori[1]))**(-1/2), (abs(autovalori[2]))**(-1/2)]
            norm_diags = normalize(diags)
            
            bpy.ops.mesh.primitive_uv_sphere_add(radius = raggioSfera * scaleFac, location=vertsAvg, scale=(norm_diags[0],norm_diags[1],norm_diags[2]))
            bpy.ops.object.shade_smooth()
            pill = bpy.context.active_object
            pill.name = obj.name + '_pill'
            assignMaterialLinear(pill, planarity, pillsColors)                                         
            blenderT_Matrix = mathutils.Matrix((eigvecs[0], eigvecs[1], eigvecs[2]))
            pill.rotation_euler = blenderT_Matrix.to_euler('XYZ')
            
            """Blender crea automaticamente gli oggetti in una master collection, ma io voglio che gli oggetti stiano nella collection che ho creato.
            Se non gestissi questa cosa, le "pillole" sarebbero collegati contemporaneamente a più collection, pertanto devo salvare tutte le vecchie collection
            alle quali l'oggetto è collegato, collegarlo alla nuova, e poi scollegarlo alle vecchie."""

            old_pill_coll = pill.users_collection
            myCol.objects.link(pill)
            for ob in old_pill_coll:
                ob.objects.unlink(pill)


#======================================================================
#                 CLASSI DI GESTIONE DELLA UI BLENDER
# (tutto il codice utilizzato nel main verrà chiamati in queste classi
#  necessarie all'aspetto grafico dell'add-on)
#======================================================================

# --- CLASSI OPERATORE ---

# Incapsula l'intero main all'interno di un operatore Blender
class MainFlowOperator(bpy.types.Operator):
    bl_idname = "scene.main_flow"
    bl_label = "Calculate!"
    bl_description = "Calculate the main geometric characteristics of the meshes, including the inertia tensor"

    def execute(self, context):
        main(context)
        return {'FINISHED'}

# Incapsula la creazione del file CSV all'interno di un operatore Blender
class CreateCSVOperator(bpy.types.Operator):
    bl_idname = "scene.select_dir"
    bl_label = "Export to CSV"
    bl_description = "Export the scanned data into a CSV file"
    filepath = bpy.props.StringProperty(subtype="DIR_PATH")
 
    def execute(self, context):
        writeCSV(self.filepath)
        return {'FINISHED'}
 
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

# --- CLASSI UI --- 

# Crea un elemento UI, in questo caso un pannello, che a sua volta
# contiene un pulsante che richiama l'operatore sopra definito
class LayoutPanel(bpy.types.Panel):
    """Creates a Panel in the scene context of the properties editor"""
    bl_label = "Point Cloud Stats"
    bl_idname = "SCENE_PT_layout"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "scene"              

    def draw(self, context):
        layout = self.layout

        scene = context.scene

        row = layout.row()
        column = layout.column()
        column.scale_y = 2.0
        row.operator("scene.main_flow")
        row.operator("scene.select_dir")

# Crea un elemento UI: visibile direttamente sulla scena 3D
# esso conterrà tutte le informazioni geometriche sulle rocce
class View3DPanel:
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "Point Cloud Stats"

    @classmethod
    def poll(cls, context):
        return (context.object is not None)

# Crea un elemento UI: visibile direttamente sulla scena 3D
# schede contenenti le informazioni geometriche
class InfoTab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_info"
    bl_label = "Info"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 7))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))
            row = self.layout.row()
            row.prop(obj, '["Object Type"]')
            meshesData[obj.name]['objType'] = obj['Object Type']

class VAvgTab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_vertex_avg"
    bl_label = "Vertices average"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 8, 11))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class VVarTab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_vertex_var"
    bl_label = "Vertices standard deviation"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 11, 14))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class NAvgTab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_normal_avg"
    bl_label = "Normals average"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 14, 17))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class NVarTab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_normal_var"
    bl_label = "Normals standard deviation"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 17, 20))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class EigValTab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_eigval"
    bl_label = "Eigenvalues"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 20, 23))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class EigVec1Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_eigvec1"
    bl_label = "Eigenvector 1"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 23, 26))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class EigVec2Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_eigvec2"
    bl_label = "Eigenvector 2"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 26, 29))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class EigVec3Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_eigvec3"
    bl_label = "Eigenvector 3"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 29, 32))
            for k, v in meshDict.items():
                self.layout.label(text = k + ": " + str(v))

class Octant1Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant1"
    bl_label = "Octant XN YN ZN"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 32, 35))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant2Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant2"
    bl_label = "Octant XN YN ZP"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 35, 38))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant3Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant3"
    bl_label = "Octant XN YP ZN"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 38, 41))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant4Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant4"
    bl_label = "Octant XN YP ZP"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 41, 44))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant5Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant5"
    bl_label = "Octant XP YN ZN"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 44, 47))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant6Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant6"
    bl_label = "Octant XP YN ZP"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 47, 50))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant7Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant7"
    bl_label = "Octant XP YP ZN"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 50, 53))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))

class Octant8Tab(View3DPanel, bpy.types.Panel):
    bl_idname = "VIEW3D_PT_octant8"
    bl_label = "Octant XP YP ZP"

    def draw(self, context):
        global meshesData
        obj = bpy.context.active_object
        if obj.name in meshesData:
            meshDict = meshesData[obj.name]
            meshDict = dict(itertools.islice(meshDict.items(), 53, 56))
            for k, v in meshDict.items():
                if k[-1] == 'x':
                    self.layout.label(text = "x: " + str(v))
                elif k[-1] == 'y':
                    self.layout.label(text = "y: " + str(v))
                elif k[-1] == 'z':
                    self.layout.label(text = "z: " + str(v))


# --- REGISTER/UNREGISTER ---
# La funzione register rende le due classi sopra definite disponibili
# alle funzionalità già presenti in Blender, così che possano essere
# integrate nel software.
# La funzione unregister fa l'esatto contrario.
def regTabs():
    bpy.utils.register_class(InfoTab)
    bpy.utils.register_class(VAvgTab)
    bpy.utils.register_class(VVarTab)
    bpy.utils.register_class(NAvgTab)
    bpy.utils.register_class(NVarTab)
    bpy.utils.register_class(EigValTab)
    bpy.utils.register_class(EigVec1Tab)
    bpy.utils.register_class(EigVec2Tab)
    bpy.utils.register_class(EigVec3Tab)
    bpy.utils.register_class(Octant1Tab)
    bpy.utils.register_class(Octant2Tab)
    bpy.utils.register_class(Octant3Tab)
    bpy.utils.register_class(Octant4Tab)
    bpy.utils.register_class(Octant5Tab)
    bpy.utils.register_class(Octant6Tab)
    bpy.utils.register_class(Octant7Tab)
    bpy.utils.register_class(Octant8Tab)

def unregTabs():
    bpy.utils.unregister_class(InfoTab)
    bpy.utils.unregister_class(VAvgTab)
    bpy.utils.unregister_class(VVarTab)
    bpy.utils.unregister_class(NAvgTab)
    bpy.utils.unregister_class(NVarTab)
    bpy.utils.unregister_class(EigValTab)
    bpy.utils.unregister_class(EigVec1Tab)
    bpy.utils.unregister_class(EigVec2Tab)
    bpy.utils.unregister_class(EigVec3Tab)
    bpy.utils.unregister_class(Octant1Tab)
    bpy.utils.unregister_class(Octant2Tab)
    bpy.utils.unregister_class(Octant3Tab)
    bpy.utils.unregister_class(Octant4Tab)
    bpy.utils.unregister_class(Octant5Tab)
    bpy.utils.unregister_class(Octant6Tab)
    bpy.utils.unregister_class(Octant7Tab)
    bpy.utils.unregister_class(Octant8Tab)

def register():
    bpy.utils.register_class(MainFlowOperator)
    bpy.utils.register_class(CreateCSVOperator)
    bpy.utils.register_class(LayoutPanel)

    regTabs()

def unregister():
    bpy.utils.unregister_class(MainFlowOperator)
    bpy.utils.unregister_class(CreateCSVOperator)
    bpy.utils.unregister_class(LayoutPanel)
    unregTabs()

if __name__ == "__main__":
    register()
