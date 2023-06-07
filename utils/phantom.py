import os.path
import random
import numpy as np
import cv2 as cv

class PhantomGenerator():
    def __init__(self,img_size, num_limit, size_limit, val_limit):
        self.img_size = img_size
        self.size_limit = size_limit
        self.val_limit = val_limit
        self.num_graphics = np.random.randint(num_limit[0], num_limit[1])
        

    def generate(self, draw_fun = None, phantomName = None, is_add_noise = False, noise_range = (-30,30)):

        # img = np.random.randint(200,255,size=[256,256]).astype(np.uint8)
        img = (np.ones(shape=self.img_size) * 255).astype(np.uint8)

        if draw_fun is None:
            draw_func = [self.draw_rectangle, self.draw_ellipse, self.draw_triangle, self.draw_trapezoidal, self.draw_parallelogram, self.draw_quadrilateral]

        for i in range(self.num_graphics):
            img = random.choice(draw_func)(img, *self.generate_para())

        if is_add_noise:
            img = (img + np.random.randint(noise_range[0], noise_range[1], size = self.img_size)).clip(0,255).astype(np.uint8)


        # if not os.path.exists(save_dir):
        #     os.makedirs(save_dir,exist_ok=True)

        return img

    def draw_ellipse(self, img, pos, shape, angle, val):

        r = (int(shape[0]/2),int(shape[1]/2))
        return cv.ellipse(img, pos[0:2], r, angle, 0, 360, val, -1)

    def draw_rectangle(self, img, pos, shape, angle, val):
        x1 = pos[0] - shape[0] / 2
        y1 = pos[1] - shape[1] / 2
        x2 = pos[0] + shape[0] / 2
        y2 = pos[1] + shape[1] / 2

        pt1 = [int(x1), int(y1)]
        pt2 = [int(x2), int(y1)]
        pt3 = [int(x2), int(y2)]
        pt4 = [int(x1), int(y2)]

        pts = [pt1, pt2, pt3, pt4]
        r_pts = []

        for pt in pts:
            x = np.int32((pt[0] - pos[0]) * np.cos(angle) + (pt[1] - pos[1]) * np.sin(angle) + pos[0])
            y = np.int32((pt[1] - pos[1]) * np.cos(angle) - (pt[0] - pos[0]) * np.sin(angle) + pos[1])
            r_pts.append([x, y])

        r_pts = np.asarray(r_pts, dtype=np.int32)

        return cv.fillPoly(img, [r_pts], val)

    def cal_distance(self,pt1,pt2):

        p1 = np.asarray(pt1)
        p2 = np.asarray(pt2)

        return np.sqrt(np.sum((p1 - p2)**2))
    def draw_triangle(self,img,pos,shape,angle,val):

        pt1 = [pos[0] + shape[0]/2,pos[1] + shape[1]/2]
        pt2 = [pos[0] - shape[2]/2,pos[1] - shape[3]/2]
        pt3 = [pos[0] + shape[4]/2,pos[1] - shape[5]/2]

        pts = [pt1, pt2, pt3]

        r_pts = []

        for pt in pts:
            x = np.int32((pt[0] - pos[0]) * np.cos(angle) + (pt[1] - pos[1]) * np.sin(angle) + pos[0])
            y = np.int32((pt[1] - pos[1]) * np.cos(angle) - (pt[0] - pos[0]) * np.sin(angle) + pos[1])
            r_pts.append([x, y])

        r_pts = np.asarray(r_pts, dtype=np.int32)

        return cv.fillPoly(img, [r_pts], val)

    def draw_parallelogram(self, img, pos, shape, angle, val):

        x1 = pos[0] - shape[0] / 2 + shape[2]
        y1 = pos[1] - shape[1] / 2
        x2 = pos[0] + shape[0] / 2 + shape[2]
        y2 = pos[1] + shape[1] / 2

        pt1 = [int(x1), int(y1)]
        pt2 = [int(x2), int(y1)]
        pt3 = [int(x2), int(y2)]
        pt4 = [int(x1), int(y2)]

        pos = pos + shape[2]

        pts = [pt1, pt2, pt3, pt4]
        r_pts = []

        for pt in pts:
            x = np.int32((pt[0] - pos[0]) * np.cos(angle) + (pt[1] - pos[1]) * np.sin(angle) + pos[0])
            y = np.int32((pt[1] - pos[1]) * np.cos(angle) - (pt[0] - pos[0]) * np.sin(angle) + pos[1])
            r_pts.append([x, y])

        r_pts = np.asarray(r_pts, dtype=np.int32)

        return cv.fillPoly(img, [r_pts], val)

    def draw_trapezoidal(self,img, pos, shape, angle, val):

        pt1 = [int(pos[0] - shape[0]/2), int(pos[1] - shape[4]/2)]
        pt2 = [int(pos[0] + shape[1]/2), int(pos[1] - shape[4]/2)]
        pt3 = [int(pos[0] + shape[2]/2), int(pos[1] + shape[4]/2)]
        pt4 = [int(pos[0] - shape[3]/2), int(pos[1] + shape[4]/2)]

        pts = [pt1, pt2, pt3, pt4]

        r_pts = []

        for pt in pts:
            x = np.int32((pt[0] - pos[0]) * np.cos(angle) + (pt[1] - pos[1]) * np.sin(angle) + pos[0])
            y = np.int32((pt[1] - pos[1]) * np.cos(angle) - (pt[0] - pos[0]) * np.sin(angle) + pos[1])
            r_pts.append([x, y])

        r_pts = np.asarray(r_pts, dtype=np.int32)

        return cv.fillPoly(img, [r_pts], val)

    def draw_quadrilateral(self,img, pos, shape, angle, val):

        pt1 = [int(pos[0] - shape[0]/2), int(pos[1] - shape[1]/2)]
        pt2 = [int(pos[0] + shape[2]/2), int(pos[1] - shape[3]/2)]
        pt3 = [int(pos[0] + shape[4]/2), int(pos[1] + shape[5]/2)]
        pt4 = [int(pos[0] - shape[6]/2), int(pos[1] + shape[7]/2)]

        pts = [pt1, pt2, pt3, pt4]

        r_pts = []

        for pt in pts:
            x = np.int32((pt[0] - pos[0]) * np.cos(angle) + (pt[1] - pos[1]) * np.sin(angle) + pos[0])
            y = np.int32((pt[1] - pos[1]) * np.cos(angle) - (pt[0] - pos[0]) * np.sin(angle) + pos[1])
            r_pts.append([x, y])

        r_pts = np.asarray(r_pts, dtype=np.int32)

        return cv.fillPoly(img, [r_pts], val)

    def generate_para(self):

        pos = tuple(np.random.randint(0, min(self.img_size[0],self.img_size[1]), size=10))
        shape = tuple(np.random.randint(self.size_limit[0], self.size_limit[1], size=10))
        angle = np.random.randint(0,360)
        val = np.random.randint(self.val_limit[0], self.val_limit[1])

        return pos, shape, angle, val


if __name__ == '__main__':
    g_img = PhantomGenerator([512,512], [1, 20], size_limit = [64, 256], val_limit = [0,128])
    g_img.generate(is_add_noise=False)
    #
    # shape = np.random.randint(10,size=10)
    #
    # print(shape)
    # print(shape[0:2]**2)

    # img = cv.imread("/home/xuepeng/ultrasound/plane_wave_beamforming/gengerate_dataset/picmus_reconstruct_img.bmp")
    # print(img.shape)
    # img = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
    # print(img.shape)
    # print(np.min(img),np.max(img))
    # print(np.random.randint(255,size = (20,20)))

    # img = np.random.randint(200,255,size=[256,256]).astype(np.uint8)
    # # img = cv.ellipse(img, (50, 50), (10, 20), 0, 0, 360, (170, 170, 170), -1)
    # # img = g_img.draw_ellipse(img,(50,50),(10,15,20),30,170)
    #
    # cv.imwrite("img.jpg",img)
    # cv.imshow('ellipse',img)
    # cv.waitKey(0)
    # shape = (1,2,3,4,5)
    # print(shape[0:2])
