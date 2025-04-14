import math
import numpy as np
from stl import mesh


# Função para criar uma base regular com a dimensão do carimbo

def createRegularBaseSolid(i_col, i_row, voxel_mm, baseT_mm, x_offset, y_offset, z_offset):
    print("Base Creator --- Calling CreateRegularBaseSolid:")
    print("> i_col",i_col)
    print("> i_row",i_row)
    print("> voxel_mm",voxel_mm)
    print("> baseT_mm",baseT_mm)
    print("> x_offset",x_offset)
    print("> y_offset",y_offset)
    print("> z_offset",z_offset)
    
    vertices = np.array([\
        [0+x_offset, 0+y_offset, 0+z_offset],
        [0+x_offset, i_row*voxel_mm+y_offset, 0+z_offset],
        [i_col*voxel_mm+x_offset, 0+y_offset, 0+z_offset],
        [i_col*voxel_mm+x_offset, i_row*voxel_mm+y_offset, 0+z_offset],
        [0+x_offset, 0+y_offset, -baseT_mm+z_offset],
        [0+x_offset, i_row*voxel_mm+y_offset, -baseT_mm+z_offset],
        [i_col*voxel_mm+x_offset, 0+y_offset, -baseT_mm+z_offset],
        [i_col*voxel_mm+x_offset, i_row*voxel_mm+y_offset, -baseT_mm+z_offset]])
    faces = np.array([\
        [0,1,2],
        [1,3,2],
        [4,1,0],
        [4,5,1],
        [3,1,7],
        [5,7,1],
        [4,7,5],
        [7,4,6],
        [3,6,2],
        [7,6,3],
        [6,0,2],
        [4,0,6]])
    solid = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            solid.vectors[i][j] = vertices[f[j],:]
    return solid

# 24/06/2021 --- Version 1 ---
# Goal: Add Reference Number to Base Assay

def createIncompleteBase(i_col, i_row, voxel_mm, baseT_mm, y_offset, z_offset):
    vertices = np.array([\
        [0, y_offset, 0+z_offset],
        [0, i_row*voxel_mm+y_offset, 0+z_offset],
        [i_col*voxel_mm, y_offset, 0+z_offset],
        [i_col*voxel_mm, i_row*voxel_mm+y_offset, 0+z_offset],
        [0, y_offset, -baseT_mm+z_offset],
        [0, i_row*voxel_mm+y_offset, -baseT_mm+z_offset],
        [i_col*voxel_mm, y_offset, -baseT_mm+z_offset],
        [i_col*voxel_mm, i_row*voxel_mm+y_offset, -baseT_mm+z_offset]])
    faces = np.array([\
        [4,1,0],
        [4,5,1],
        [3,1,7],
        [5,7,1],
        [4,7,5],
        [7,4,6],
        [3,6,2],
        [7,6,3],
        [6,0,2],
        [4,0,6]])
    cube = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            cube.vectors[i][j] = vertices[f[j],:]
    return cube


# Main Function
def createBaseTopWithReferenceHoles(reference, b_voxel, letter_deep, i_col, i_row, l_voxel, pos_x, pos_y, z_offset, axis):
    
    global total_mesh_vertices, total_mesh_faces, index, b_voxel_mm, l_voxel_mm
    total_mesh_vertices = []
    total_mesh_faces = []
    index = 0
    b_voxel_mm = b_voxel
    l_voxel_mm = l_voxel
    
    topright_y_px = 0
    bottomleft_x_px = 0
    bottomleft_y_px = 0
    topright_x_px = 0
    if axis=="AxisX":
        # AXIS X
        topright_y_px = pos_y
        bottomleft_x_px = pos_x
        bottomleft_y_px = topright_y_px+7
        topright_x_px = bottomleft_x_px+len(reference)*7
    else:
        # AXIS Y
        topright_y_px = pos_y
        bottomleft_x_px = pos_x
        bottomleft_y_px = topright_y_px+len(reference)*7
        topright_x_px = bottomleft_x_px+7 
    
    # Side L
    for row_y_px in range(topright_y_px,bottomleft_y_px):
        total_mesh_vertices.append([bottomleft_x_px*l_voxel_mm, row_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([bottomleft_x_px*l_voxel_mm, (row_y_px+1)*l_voxel_mm, z_offset])
        total_mesh_vertices.append([0,row_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([0,(row_y_px+1)*l_voxel_mm, z_offset])
        total_mesh_faces.append([index+2,index+1,index+0])
        total_mesh_faces.append([index+3,index+1,index+2])
        index += 4
    # Side R
    for row_y_px in range(topright_y_px,bottomleft_y_px):
        total_mesh_vertices.append([topright_x_px*l_voxel_mm, row_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([topright_x_px*l_voxel_mm, (row_y_px+1)*l_voxel_mm, z_offset])
        total_mesh_vertices.append([i_col*b_voxel_mm,row_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([i_col*b_voxel_mm,(row_y_px+1)*l_voxel_mm, z_offset])
        total_mesh_faces.append([index+0,index+1,index+2])
        total_mesh_faces.append([index+2,index+1,index+3])
        index += 4
    # Side B
    for col_x_px in range(bottomleft_x_px,topright_x_px):
        total_mesh_vertices.append([col_x_px*l_voxel_mm, topright_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([(col_x_px+1)*l_voxel_mm, topright_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([col_x_px*l_voxel_mm,0, z_offset])
        total_mesh_vertices.append([(col_x_px+1)*l_voxel_mm,0, z_offset])
        total_mesh_faces.append([index+0,index+1,index+2])
        total_mesh_faces.append([index+2,index+1,index+3])
        index += 4
    # Side T
    for col_x_px in range(bottomleft_x_px,topright_x_px):
        total_mesh_vertices.append([col_x_px*l_voxel_mm, bottomleft_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([(col_x_px+1)*l_voxel_mm, bottomleft_y_px*l_voxel_mm, z_offset])
        total_mesh_vertices.append([col_x_px*l_voxel_mm,i_row*b_voxel_mm, z_offset])
        total_mesh_vertices.append([(col_x_px+1)*l_voxel_mm,i_row*b_voxel_mm, z_offset])
        total_mesh_faces.append([index+2,index+1,index+0])
        total_mesh_faces.append([index+3,index+1,index+2])
        index += 4
    # Big Square 1
    total_mesh_vertices.append([bottomleft_x_px*l_voxel_mm, bottomleft_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([0, bottomleft_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([0, i_row*b_voxel_mm, z_offset])
    total_mesh_vertices.append([bottomleft_x_px*l_voxel_mm, i_row*b_voxel_mm, z_offset])
    total_mesh_faces.append([index+0,index+1,index+2])
    total_mesh_faces.append([index+2,index+3,index+0])
    index += 4
    # Big Square 2
    total_mesh_vertices.append([bottomleft_x_px*l_voxel_mm, topright_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([0, topright_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([0, 0, z_offset])
    total_mesh_vertices.append([bottomleft_x_px*l_voxel_mm, 0, z_offset])
    total_mesh_faces.append([index+2,index+1,index+0])
    total_mesh_faces.append([index+0,index+3,index+2])
    index += 4
    # Big Square 3
    total_mesh_vertices.append([topright_x_px*l_voxel_mm, topright_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([i_col*b_voxel_mm, topright_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([i_col*b_voxel_mm, 0, z_offset])
    total_mesh_vertices.append([topright_x_px*l_voxel_mm, 0, z_offset])
    total_mesh_faces.append([index+0,index+1,index+2])
    total_mesh_faces.append([index+2,index+3,index+0])
    index += 4
    # Big Square 4
    total_mesh_vertices.append([topright_x_px*l_voxel_mm, bottomleft_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([i_col*b_voxel_mm, bottomleft_y_px*l_voxel_mm, z_offset])
    total_mesh_vertices.append([i_col*b_voxel_mm, i_row*b_voxel_mm, z_offset])
    total_mesh_vertices.append([topright_x_px*l_voxel_mm, i_row*b_voxel_mm, z_offset])
    total_mesh_faces.append([index+2,index+1,index+0])
    total_mesh_faces.append([index+0,index+3,index+2])
    index += 4
    # Generate Stamp Reference
    start_x_px = bottomleft_x_px
    start_y_px = topright_y_px
    start_z_px = z_offset
    for letter in reference:
        if letter == '1':
            drawOne(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '2':
            drawTwo(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '3':
            drawThree(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '4':
            drawFour(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '5':
            drawFive(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '6':
            drawSix(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '7':
            drawSeven(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '8':
            drawEight(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '9':
            drawNine(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '0':
            drawZero(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == '.':
            drawPoint(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'a':
            drawLetterA(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'b':
            drawLetterB(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'c':
            drawLetterC(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'd':
            drawLetterD(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'e':
            drawLetterE(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'f':
            drawLetterF(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'g':
            drawLetterG(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'h':
            drawLetterH(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'i':
            drawLetterI(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'j':
            drawLetterJ(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'k':
            drawLetterK(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'l':
            drawLetterL(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'm':
            drawLetterM(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'n':
            drawLetterN(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'o':
            drawLetterO(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'p':
            drawLetterP(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'q':
            drawLetterQ(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'r':
            drawLetterR(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 's':
            drawLetterS(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 't':
            drawLetterT(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'u':
            drawLetterU(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'v':
            drawLetterV(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'w':
            drawLetterW(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'x':
            drawLetterX(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'y':
            drawLetterY(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        elif letter == 'z':
            drawLetterZ(letter_deep, start_x_px, start_y_px, start_z_px, axis)
        if axis=="AxisX":
            start_x_px += 7
        else:
            start_y_px += 7
        
    # Mesh
    vertices = np.array(total_mesh_vertices)
    faces = np.array(total_mesh_faces)
    volume = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            volume.vectors[i][j] = vertices[f[j],:]
    return volume


# POINT
def drawPoint(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0],
                     [0,0,1,1,1,0,0],
                     [0,0,1,1,1,0,0],
                     [0,0,1,1,1,0,0],
                     [0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 0 ZERO
def drawZero(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 1 ONE
def drawOne(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,1,1,0,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 2 TWO
def drawTwo(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 3 THREE
def drawThree(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,1,0],
                     [0,0,1,1,1,1,0],
                     [0,0,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 4 FOUR
def drawFour(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,0,0,1,0],
                     [0,0,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 5 FIVE
def drawFive(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 6 SIX
def drawSix(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 7 SEVEN
def drawSeven(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,0,0,0,0],
                     [0,0,1,1,0,0,0],
                     [0,0,0,1,1,0,0],
                     [0,0,0,0,1,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 8 EIGHT
def drawEight(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# 9 NINE
def drawNine(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,0,0,1,0],
                     [0,0,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# A LETTER
def drawLetterA(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# B LETTER
def drawLetterB(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,1,0,0,0],
                     [0,1,1,1,0,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# C LETTER
def drawLetterC(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# D LETTER
def drawLetterD(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,0,0],
                     [0,1,0,0,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,1,1,0],
                     [0,1,1,1,1,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# E LETTER
def drawLetterE(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# F LETTER
def drawLetterF(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# G LETTER
def drawLetterG(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,0,0],
                     [0,1,0,0,1,0,0],
                     [0,1,0,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,1,1,1,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# H LETTER
def drawLetterH(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# I LETTER
def drawLetterI(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,1,1,1,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,1,1,1,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# J LETTER
def drawLetterJ(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,0,1,0],
                     [0,0,0,0,0,1,0],
                     [0,0,0,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# K LETTER
def drawLetterK(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,1,1,0],
                     [0,1,0,1,1,0,0],
                     [0,1,1,1,0,0,0],
                     [0,1,0,1,1,0,0],
                     [0,1,0,0,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# L LETTER
def drawLetterL(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,1,0,0,0,0,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# M LETTER
def drawLetterM(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,1,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,1,0,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# N LETTER
def drawLetterN(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,1,1,0],
                     [0,1,0,1,1,1,0],
                     [0,1,0,1,0,1,0],
                     [0,1,1,1,0,1,0],
                     [0,1,1,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# O LETTER
def drawLetterO(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,0,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,0,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# P LETTER
def drawLetterP(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,0,0,0,0],
                     [0,1,1,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,1,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# Q LETTER
def drawLetterQ(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,0,0,0,0],
                     [0,0,1,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# R LETTER
def drawLetterR(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,1,1,0],
                     [0,1,1,1,1,0,0],
                     [0,1,0,0,1,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# S LETTER
def drawLetterS(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,1,1,0],
                     [0,1,1,1,1,1,0],
                     [0,1,1,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# T LETTER
def drawLetterT(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# U LETTER
def drawLetterU(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# V LETTER
def drawLetterV(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,1,1,1,0,0],
                     [0,1,1,0,1,1,0],
                     [0,1,1,0,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# W LETTER
def drawLetterW(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,1,0,1,0,0],
                     [0,0,1,1,1,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,0,1,0,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# X LETTER
def drawLetterX(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,0,0,0,1,0],
                     [0,1,1,0,1,1,0],
                     [0,0,1,1,1,0,0],
                     [0,1,1,0,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# Y LETTER
def drawLetterY(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,0,1,0,0,0],
                     [0,0,1,1,1,0,0],
                     [0,1,1,0,1,1,0],
                     [0,1,0,0,0,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# Z LETTER
def drawLetterZ(drop, start_x_px, start_y_px, start_z_px, rotation_axis):
    
    num = np.matrix([[0,0,0,0,0,0,0],
                     [0,1,1,1,1,1,0],
                     [0,1,1,0,0,0,0],
                     [0,0,1,1,1,0,0],
                     [0,0,0,0,1,1,0],
                     [0,1,1,1,1,1,0],
                     [0,0,0,0,0,0,0]])
    if rotation_axis=="AxisX":
        loopMatrix2Draw(num, drop, start_x_px, start_y_px, start_z_px)
    else:
        loopMatrix2Draw(np.rot90(num,k=3), drop, start_x_px, start_y_px, start_z_px)

# MATRIX LOOP DRAWER
def loopMatrix2Draw(matrix, drop, start_x_px, start_y_px, start_z_px):
    global index
    # Loop Throw Number Matrix (Top and Bottom)
    for i in range(0,matrix.shape[0]):
        for j in range(0,matrix.shape[1]):
            if matrix[i,j]==1:
                z_coord = -drop
                total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord])
                total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord])
                total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord])
                total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord])
                total_mesh_faces.append([index+2,index+1,index+0])
                total_mesh_faces.append([index+3,index+1,index+2])
                index += 4
            else:
                z_coord = start_z_px
                total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord])
                total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord])
                total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord])
                total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord])
                total_mesh_faces.append([index+2,index+1,index+0])
                total_mesh_faces.append([index+3,index+1,index+2])
                index += 4
    # Loop Throw Number Matrix (Sides)
    z_coord_bottom = -drop
    z_coord_top = start_z_px
    for i in range(0,matrix.shape[0]):
        for j in range(0,matrix.shape[1]):
            if i+1<=matrix.shape[0]:
                if matrix[i,j]==1 and matrix[i+1,j]==0:
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_bottom])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_bottom])
                    total_mesh_faces.append([index+0,index+1,index+2])
                    total_mesh_faces.append([index+2,index+1,index+3])
                    index += 4
            if i-1>=0:
                if matrix[i,j]==1 and matrix[i-1,j]==0:
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_bottom])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_bottom])
                    total_mesh_faces.append([index+2,index+1,index+0])
                    total_mesh_faces.append([index+3,index+1,index+2])
                    index += 4 
            if j-1>=0:
                if matrix[i,j]==1 and matrix[i,j-1]==0:
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_bottom])
                    total_mesh_vertices.append([(start_x_px+j)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_bottom])
                    total_mesh_faces.append([index+2,index+1,index+0])
                    total_mesh_faces.append([index+3,index+1,index+2])
                    index += 4
            if j+1<=matrix.shape[1]:
                if matrix[i,j]==1 and matrix[i,j+1]==0:
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_top])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i+1)*l_voxel_mm, z_coord_bottom])
                    total_mesh_vertices.append([(start_x_px+j+1)*l_voxel_mm, (start_y_px+i)*l_voxel_mm, z_coord_bottom])
                    total_mesh_faces.append([index+0,index+1,index+2])
                    total_mesh_faces.append([index+2,index+1,index+3])
                    index += 4


# Public Interface -------------------------------------------------------------------                   
def generateBaseWithReference(i_col, i_row, stamp_voxel_xy, base_thickness, base_z_offset, letter_pos_x, letter_pos_y, letter_voxel_xy, letter_engraving_z, letter_reference, letter_axis):
    img_col = i_col
    img_row = i_row
    stamp_voxel_xy_mm = stamp_voxel_xy
    base_thickness_mm = base_thickness
    base_z_offset_mm = base_z_offset
    
    letter_pos_x_px = letter_pos_x
    letter_pos_y_px = letter_pos_y
    letter_voxel_xy_mm = letter_voxel_xy
    letter_engraving_z_mm = letter_engraving_z
    letter_reference_str = letter_reference
    
    base3D = createIncompleteBase(img_col,img_row,stamp_voxel_xy_mm,base_thickness_mm,0.0,base_z_offset_mm)
    
  
    reference3D = createBaseTopWithReferenceHoles(letter_reference_str,
                                                  stamp_voxel_xy_mm,
                                                  letter_engraving_z_mm,
                                                  img_col,
                                                  img_row,
                                                  letter_voxel_xy_mm,
                                                  letter_pos_x_px,
                                                  letter_pos_y_px,
                                                  base_z_offset_mm,
                                                  letter_axis)
    
    base_final_mesh = mesh.Mesh(np.concatenate([
        reference3D.data.copy(),
        base3D.data.copy(),
        ]))
    
    return base_final_mesh



# Development Test -------------------------------------------------------------------
# baseNew = generateBaseWithReference(1228,
#                                     2201,
#                                     0.02,
#                                     1.0,
#                                     0.05,
#                                     2,
#                                     2,
#                                     0.4,
#                                     0.5,
#                                     '1234.0001.000',
#                                     'AxisY')
#     
# baseNew.save('stamp_base_v3.stl')