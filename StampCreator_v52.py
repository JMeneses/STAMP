from skimage import morphology
from skimage import measure
import math
import numpy as np
from PIL import Image
from stl import mesh
import BaseCreator_v1 as bc

# Parameters:
# stampZ_mm > Upper part of the stamp height ZZ (top where will be the objective image)
# stampY_mm > Bottom part of the stamp length YY
# stampX_mm > Bottom part of the stamp width XX
# stampBorder_mm > Upper part of the stamp border, margin to each side, lateral growth (image dilation) 
# stampBaseThickness_mm > Bottom part of the stamp height
# sliceSize_mm > define layers thickness, vertical growth
# voxelSize_mm > define XY voxel size

# Variables
progress_percent = 0
stampZ_mm = 2.0 
stampY_mm = 49
stampX_mm = 74
stampBorder_mm = stampZ_mm/np.tan(60*np.pi/180) #Introducing an angle consideration, 60 degrees
stampBaseThickness_mm = 1.0
sliceSize_mm = 0.05
voxelSize_mm = 0.05

print("Initial construction parameters:")
print("stampZ_mm > ",stampZ_mm," mm")
print("stampY_mm > ",stampY_mm," mm")
print("stampX_mm > ",stampX_mm," mm")
print("stampBorder_mm (dilation) > ",stampBorder_mm," mm")
print("sliceSize_mm > ",sliceSize_mm," mm")
print("voxelSize_mm > ",voxelSize_mm," mm")
print("")
print("")

# 0. Load Image (Force Greyscale)
print("0. Welcome to STAMP (v52):")
img = Image.open('C:/Users/momen/Desktop/STAMP OpenSource/Ficheiros_Teste/RIOJA-DIAM_BW_flipped_size005mm.png',mode='r')
img = img.convert("L")
imgArray = np.array(img)
process_percent = 10
print("0.1 Loaded img size (x,y) > ",img.size," (required resolution to fullfill total size criteria x:",stampX_mm/voxelSize_mm,"y:",stampY_mm/voxelSize_mm,")")
if (img.size[0]!=stampX_mm/voxelSize_mm or img.size[1]!=stampY_mm/voxelSize_mm):
    print("ERROR >> Invalid image size, image resultion must be exactly equal to x:",stampX_mm/voxelSize_mm,"y:",stampY_mm/voxelSize_mm)
print("")
print("")

# 1. Read Image Data
img_col, img_row = imgArray.shape
print("1. Generate Ground Map (0 Outer and 11 Inner).")
img_borders = np.zeros((img_col, img_row))
for x in range(0, img_col):
    for y in range(0, img_row):
        if (imgArray[x,y]<=122):
            img_borders[x,y] = 0
        else:
            img_borders[x,y] = 11
process_percent = 20
print("")
print("")
     
# 2. Slice Union to 3D Volume (Quality Control Rules)
print("2. Slice Union to 3D Volume")





print("2.1. RULE For the selected base and voxel sizes, the recommended image resolution is: ",stampX_mm/voxelSize_mm,"x",stampY_mm/voxelSize_mm,"y")
print("2.2. RULE The image must have just one color channel (GREY) and must have been binarized (0,1) or (0,255), black means flat background, white means the top of the stamp.")
print("2.3. RULE The image must have a border padding with ",round(stampBorder_mm/voxelSize_mm)," pixels, they must have only zero value elements.")

if (stampY_mm/img_col == stampX_mm/img_row):
    # 1st Quality Criteria
    if (stampY_mm/img_row<=voxelSize_mm):
        # 2nd Quality Criteria
        voxelSize_mm = stampX_mm/img_row
    else:
        # Warning (2) - Voxel quality is inferior to the recomended quality
        print("Warning (2) - Voxel quality is inferior to the recomended quality.")
        voxelSize_mm = stampX_mm/img_row
else:
    # Error (1) - Image dimensions don't match the desired stamp ratio, advise new image resolution, can't proceed
    print("Error (1) - Image dimensions don't match the desired stamp ratio, advise new image resolution, can't proceed.",(stampY_mm/img_col),"==",(stampX_mm/img_row))





# 2.2. Calulate Dependent Variables 
sliceN = round(stampZ_mm/sliceSize_mm)
borderN = math.ceil(stampBorder_mm/voxelSize_mm)
stepN = math.ceil(borderN/sliceN)
dilationRange = range(borderN, 0, -stepN)

print("sliceN > ",sliceN," [ad]") #Number of slices in the Z direction
print("borderN > ",borderN," [ad]") #Number of radial expansions to complete the defined stampBorder_mm
print("stepN > ",stepN," [ad]") #
print("dilationRange > ",dilationRange," [ad]") #dilation by layer
print("")

# 2.3. Image Dilation (adding slices from bottom to top)
img_volume = np.zeros((img_col, img_row, sliceN))
for n in range(0,sliceN):
    if (n==0 or n==sliceN-1):
        img_volume[:,:,n] = np.zeros((img_col, img_row)) # the first and the last layer need to be empty space
    else:
        if (n==sliceN-2):
            img_volume[:,:,n] = img_borders # the objective stamp image
        else:
            if (n>=len(dilationRange)):
                img_volume[:,:,n] = morphology.dilation(img_borders,np.ones((dilationRange[-1],dilationRange[-1]))) # adding layer without dilation
            else:
                img_volume[:,:,n] = morphology.dilation(img_borders,np.ones((dilationRange[n],dilationRange[n]))) # adding layer with dilation





# sliceN min 3 to be valid




process_percent = 30




# Marching cubes
print("3. Marching Cubes Algorithm")
verts, faces, normals, values = measure.marching_cubes(img_volume) #removi o segundo parâmetro 10

# new size
print("4. Size Conversion")
for v in verts:
    v[0] = v[0]*voxelSize_mm
    v[1] = v[1]*voxelSize_mm
    v[2] = v[2]*sliceSize_mm

# Numpy - STL
print("5. STL Model Mesh Creation")
model3D = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        model3D.vectors[i][j] = verts[f[j],:]
        
# Add Base Mesh To Other Mesh
print("6. STL Base Mesh Creation")
base3D = bc.createRegularBaseSolid(img_col,img_row,voxelSize_mm,stampBaseThickness_mm,0.0,0.0,sliceSize_mm)
final_mesh = mesh.Mesh(np.concatenate([
    model3D.data.copy(),
    base3D.data.copy(),
    ]))

# Write the mesh to STL file 
print("7. Save STL Model")
final_mesh.save('last_generated_3D_model_stamp_V52.stl')

print("8. The End")

