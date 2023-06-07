import os
import re
import cv2
import time
import logging
import matlab
import matlab.engine
import multiprocessing

from multiprocessing import Process
from utils.phantom import PhantomGenerator

mateng = matlab.engine.start_matlab('-nodesktop -nodisplay')
mateng.cd(r"C:\Users\Administrator\Documents\MATLAB\CPWC\matlab_script")
# mateng.addpath(r'C:\Users\hrzy\Documents\MATLAB\CPWC\matlab_script\picmus')
# mateng.addpath(r'C:\Users\hrzy\Documents\MATLAB\CPWC\matlab_script\Field_II_ver_3_30_linux')
# mateng.addpath(r'C:\Users\hrzy\Documents\MATLAB\CPWC\matlab_script\ustb')


def mat_reconstruct_img(rawDataDir, usImageDir, phantomID):

    raw_data_path = os.path.join(rawDataDir, f'rfdata_{phantomID}')
    usimage_path = os.path.join(usImageDir, f'usimage_{phantomID}.png')

    mateng.reconstruct_img(raw_data_path, usimage_path)


def generate_phantom_img(phantomImgDir, phantomID):

    phantomName = os.path.join(phantomImgDir, f'phantom_image_{phantomID}.png')
    g_img = PhantomGenerator([512,512], [1, 24], size_limit = [64, 256], val_limit = [32,128])
    phantomImg = g_img.generate(phantomName=phantomName, is_add_noise=False)

    cv2.imwrite(phantomName, phantomImg)

def mat_generate_phantom(phantomImgDir, phantomMatDir, phantomID, N):

    phantom_img_path = os.path.join(phantomImgDir, f'phantom_image_{phantomID}.png')
    phantom_matrix_path = os.path.join(phantomMatDir, f'phantom_matrix_{phantomID}.hdf5')

    mateng.make_scatters(phantom_img_path, phantom_matrix_path, N)


def mat_sim(processid, phantomMatDir, rawDataDir, phantomID, start, end):
    start_angle = float(start)
    end_angle = float(end)
    print(f'sub process: {processid}({os.getpid()}) start...')

    phantom_matrix_path = os.path.join(phantomMatDir, f'phantom_matrix_{phantomID}.hdf5')
    save_dir = os.path.join(rawDataDir, f'rfdata_{phantomID}')

    time.sleep(10)

    mateng.sim_cpwc(start_angle, end_angle, phantom_matrix_path, save_dir)


def fork_process(phantomMatDir, rawDataDir, phantomID):
    try:
        start_time = time.time()
        iProcesses = []

        for i in range(1, 16):

            end_angle = 5 * i
            start_angle = end_angle - 4
            processid = 'SUB PROCESS #{:>2}'.format(i)
            iprocess = Process(target = mat_sim, args=(processid, phantomMatDir, rawDataDir, phantomID, start_angle, end_angle))
            iProcesses.append(iprocess)

        for iprocess in iProcesses:
            iprocess.start()
            
        for iprocess in iProcesses:
            iprocess.join()

        end_time = time.time()
        print(f'consumed time: {end_time - start_time}')

    except KeyboardInterrupt:
        print('KeyboardInterrupt....')
        # mateng.quit()
        for iprocess in iProcesses:
            iprocess.terminate()
    # else:
    #     for iprocess in iProcesses:
    #         iprocess.terminate()
    # finally:
    #     # mateng.quit()
    #     for iprocess in iProcesses:
    #         iprocess.terminate()


def config(root):
    phantomImgDir = 'phantom_image'
    phantomMatDir = 'phantom_matrix'
    rawDataDir = 'rfdata'
    usImageDir = 'usimage'

    phantomImgDir = os.path.join(root, phantomImgDir)
    phantomMatDir = os.path.join(root, phantomMatDir)
    rawDataDir = os.path.join(root, rawDataDir)
    usImageDir = os.path.join(root, usImageDir)

    if not os.path.exists(phantomImgDir):
        os.mkdir(phantomImgDir)

    if not os.path.exists(phantomMatDir):
        os.mkdir(phantomMatDir)

    if not os.path.exists(rawDataDir):
        os.mkdir(rawDataDir)

    if not os.path.exists(usImageDir):
        os.mkdir(usImageDir)

    return phantomImgDir, phantomMatDir, rawDataDir, usImageDir


def findMaxID(fileName):
    items = os.listdir(fileName)

    maxID = 0
    for item in items:
        id = int(re.findall('\d+', item)[0])
        if id > maxID:
            maxID = id
    
    return maxID



if __name__ == '__main__':
    # name = 'CPWC20230525'

    fileName = r'C:\Users\Administrator\Documents\MATLAB\CPWC\dataset\train'
    N = 1e5

    phantomImgDir, phantomMatDir, rawDataDir, usImageDir = config(fileName)

    phantomID = 0
    
    print(f'Start simulate phantom from ID-{phantomID}')

    multiprocessing.set_start_method('spawn')
    print(f'master process ({os.getpid()}) start...')

    PHANTOM_NUM = 0
    epoch = 0

    print('==>generate phantom original image...')
    generate_phantom_img(phantomImgDir, phantomID)

    print('==>generate phantom matrix...')
    mat_generate_phantom(phantomImgDir, phantomMatDir, phantomID, N)

    print('==>simulate cpwc with multiprocessing...')
    fork_process(phantomMatDir, rawDataDir, phantomID)

    print('==>cpwc simulating finished...')
    
            
    print('\n==>reconstruct ultrasound simulated imgae...')
    mat_reconstruct_img(rawDataDir, usImageDir, phantomID)


    mateng.quit()