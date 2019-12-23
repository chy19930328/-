import numpy as np
import cv2
import math
import bayes
def create_6_box(im,center,the_long,the_color_choose_index):
    point1 = [int(center[1] - math.sqrt(3)*the_long/2),int(center[0] - the_long/2)]
    point2 = [int(center[1] - math.sqrt(3)*the_long/2),int(center[0] + the_long/2)]
    point3 = [center[1],center[0]+the_long]
    point4=  [int(center[1] + math.sqrt(3)*the_long/2),int(center[0] + the_long/2)]
    point5 = [int(center[1] + math.sqrt(3)*the_long/2),int(center[0] - the_long/2)]
    point6 = [int(center[1]),int(center[0] -the_long)]
    a = np.array([[point1,point2,point3,point4,point5,point6]],dtype=np.int32)
    a =a[:,:,(1,0)]
    the_color =np.array( [[0,255,255],[0,0,255],[255,0,0],[255,255,0],[255,255,255]]).astype(np.int)
    the_choose_c =np.array( (the_color[the_color_choose_index][0],the_color[the_color_choose_index][1],the_color[the_color_choose_index][2]))
    color = tuple([int(x) for x in the_choose_c]) 
    cv2.fillPoly(im, a, color=color)
    cv2.polylines(im, a,True,(255,255,255))
    return im

message =np.load("Tdata.npy").astype(np.int)
the_long = 12
the_first_center = [int(math.sqrt(3)*the_long/2),the_long]
the_x_jiange = int(the_long*3/2)
the_y_jiange = int(math.sqrt(3)*the_long)
h,w = message.shape

the_image_height = int(the_y_jiange*message.shape[0]+math.sqrt(3)*the_long)+1
the_image_width = int(the_x_jiange*message.shape[1]+2*the_long)+1

image = np.zeros((the_image_height,the_image_width,3)).astype(np.uint8)

for i in range(h):
    for j in range(w):
        the_init_jiange_x = the_first_center[1]
        the_init_jiange_y = the_first_center[0]
        if j%2==1:
            the_init_jiange_y+=the_init_jiange_y
#            continue
        
        the_x = the_init_jiange_x+ the_x_jiange*j
        the_y = the_init_jiange_y+the_y_jiange*i
        the_color_index = message[i,j]

        create_6_box(image,[the_x,the_y],the_long,the_color_index)
#image = cv2.resize(image,(image.shape[1]/4),image.shape[0]/4))
#cv2.imshow("dfa",image)
#cv2.waitKey(0)


the_path =bayes.Path
for i in range(len(the_path)-1):
    the_one = [int(the_path[i].split("+")[0]),int(the_path[i].split("+")[1])]
    the_two = [int(the_path[i+1].split("+")[0]),int(the_path[i+1].split("+")[1])]
    
    kki = the_one[0]
    kkj =the_one[1]
    the_init_jiange_x = the_first_center[1]
    the_init_jiange_y = the_first_center[0]
    if kkj%2==1:
        the_init_jiange_y+=the_init_jiange_y
    
    the_x1 = the_init_jiange_x+ the_x_jiange*kkj
    the_y1 = the_init_jiange_y+the_y_jiange*kki

    
    kki = the_two[0]
    kkj =the_two[1]
    the_init_jiange_x = the_first_center[1]
    the_init_jiange_y = the_first_center[0]
    if kkj%2==1:
        the_init_jiange_y+=the_init_jiange_y
    
    the_x2 = the_init_jiange_x+ the_x_jiange*kkj
    the_y2 = the_init_jiange_y+the_y_jiange*kki
    lineType = 4

    image = cv2.arrowedLine(image, (the_x1,the_y1), (the_x2,the_y2), (4,100,4),2,1,0,0.3)
#    cv2.line(image, (the_x1,the_y1), (the_x2,the_y2), (255,255,0), 2, lineType)

cv2.imshow("dfa",image)
cv2.waitKey(0)