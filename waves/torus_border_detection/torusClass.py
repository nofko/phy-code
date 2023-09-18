import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import medfilt,argrelextrema
import os
from os.path import exists
import sys

class torus:
    """Class used for analyzing the torus, contains all neccessary functions and variables

    Input:
        filenames   - list of 3 elements, file locations [video,calibration_full, calibration_empty]
        r_phys      - float, physical radius of the interface we are analyzing
        r_in        - int, radius of inner crop circle
        r_out       - int, radius of outer crop circle
    """

    def __init__(self,filename,r_phys,r_in,r_out):

        self.videoLoc=filename

        self.showVideo=1
        
        self.r_phys=r_phys
        self.r_in=r_in
        self.r_out=r_out

        self.seuil=5
        self.kerSize=7

        self.xCenter=0
        self.yCenter=0

        self.D1=1850
        self.D2=5
        self.kA=5
        self.thetaLen=3801
        self.thetaLin=np.linspace(0,2*np.pi,num=self.thetaLen,endpoint=True)

        self.g_eff=9.81*np.sin(4.5*np.pi/180)
        self.sigma_eff=0.045
        self.c_eff=0.07
        self.rho=1000

        self.cmin=7.4
        self.cmax=15

        self.calibImg=0
        self.closeCont=0
        self.percInit=0.9

        return

    def openVideo(self):
        """Opens the video and calibration images, and stores number of frames and fps"""

        self.video=cv2.VideoCapture(self.videoLoc,0)

        #self.calibFull=cv2.imread(self.calibFullLoc,0)

        #self.calibEmpty=cv2.imread(self.calibEmptyLoc,0)

        self.numFrames=int(self.video.get(cv2.CAP_PROP_FRAME_COUNT))
        self.fps=int(self.video.get(cv2.CAP_PROP_FPS))

        return

    def runVideo(self):
        """Loops over the video to give a final output signal

        Input:
            video - the video which we want to analyze

        Output:
            numpy array - (time,theta) numpy array containng the signal in time and space

        """

        if self.calibImg==1:
            self.calibrate(self.calibFull)
            #plt.imshow(self.calibFull)
        else:
            ret,frame=self.video.read()
            self.calibrate(frame)
            #plt.imshow(cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY))

        #self.seuil=int(input("Seuil: "))
            
        self.signal=np.zeros((self.numFrames,self.thetaLen))
        self.rawSig=np.zeros((self.numFrames,self.thetaLen))

        
        i=0
        while(self.video.isOpened()):

            ret, frame = self.video.read()

            if ret == True:

                border=self.getBorder(frame,0)
                tempSig=self.getSignal(border)
                self.signal[i,:]=tempSig-self.calibSig

                self.rawSig[i,:]=tempSig

                i+=1

                #os.system('clear')
                sys.stdout.flush()
                sys.stdout.write("Progress: "+str(int(100*i/self.numFrames))+" %\r")
            else:
                break

        print("Done")
        self.video.release()

        return

    def getBorder(self,image,calib):
        """Find the x and y coordinates, around the center of the image, of the circular contour

        Input:
            image - an image, or a frame from a video
            calib - integer, 1 if the image is a calibration image

        Output:
            numpy array - 1st element = list, x coordinates
                        - 2nd element = list, y coordinates

        Key parameters:
            r_inner - defined globally, inside crop
            r_outer - defined globally, outside crop
            seuil   - to be adjusted depending on image

        """

        height, width = image.shape[:2]

        if(len(image.shape)<3):
            gray=image
        else:
            gray=cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

        mask_inner=np.ones_like(gray)
        mask_outer=np.zeros_like(gray)

        cv2.circle(mask_inner,(int(width/2),int(height/2)),self.r_in,0,-1,8,0)
        cv2.circle(mask_outer,(int(width/2),int(height/2)),self.r_out,1,-1,8,0)

        circ=mask_inner*(1-gray)
        circ*=mask_outer

        circ=np.where(circ>self.seuil,255,0)
        circ=circ.astype(np.uint8)

        rot_rectangle = ((870, 1300), (80, 260), -11)

        box = cv2.boxPoints(rot_rectangle) 
        box = np.int0(box)

        cv2.fillPoly(circ, pts =[box], color=(255,255,255))
        
        if(self.closeCont==1):
            kernel=np.ones((self.kerSize,self.kerSize),np.uint8)
            circ=cv2.morphologyEx(circ,cv2.MORPH_CLOSE,kernel)

        contours, hierarchy = cv2.findContours(circ, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

        cMax=max(contours,key=cv2.contourArea)

        if calib==1:

            self.rect=cv2.fitEllipse(cMax)
            self.xCenter=self.rect[0][0]
            self.yCenter=self.rect[0][1]

        if self.showVideo==1:
            
            blank=np.ones((height,width,3),np.uint8)*255
            blank=cv2.cvtColor(gray,cv2.COLOR_GRAY2RGB)
            cv2.drawContours(blank, [cMax], 0, (255,0,0),3)
            
            
            #cv2.drawContours(blank, contours, -1, (255,0,0),3)

            cv2.circle(blank,(int(self.xCenter),int(self.yCenter)),int(self.r_in),(0,0,255))
            cv2.circle(blank,(int(self.xCenter),int(self.yCenter)),int(self.r_out),(0,0,255))

            #rectangle = cv2.drawContours(blank,[box],0,(0,0,255),2)
            
            cv2.namedWindow('shown',cv2.WINDOW_NORMAL)
            cv2.imshow('shown',blank)
            cv2.resizeWindow('shown',800,600)
            cv2.waitKey(2)

        shape=np.vstack(cMax).squeeze()
        xSig=shape[:,0]-self.xCenter
        ySig=shape[:,1]-self.yCenter

        return np.array([xSig,ySig])

    def getSignal(self,cont):
        """Takes in coordinates of a contour and returns the interpolated signal

        Input:
            contour - numpy array [X,Y] where X and Y are arrays of coordinates

        Output:
            [theta,radius] - two elements, each is a list of theta (radius) coordinates
        """
        xar=cont[0].astype(float)
        yar=cont[1].astype(float)
        car=np.divide(yar,xar,np.ones_like(xar)*np.pi,where=xar!=0)

        radialPos=np.sqrt(xar**2+yar**2)
        thetaPos=np.arctan(car)
        #thetaPos=np.where(xar==0.,np.pi,np.arctan(yar/xar))

        #for i in range(len(cont[0])):
        #    if cont[0][i]<0:
        #        thetaPos[i]+=np.pi
        #    if cont[0][i]>0 and cont[1][i]<0:
        #        thetaPos[i]+=2*np.pi

        thetaPos=np.where(xar<0,thetaPos+np.pi,thetaPos)
        thetaPos=np.where((xar>0) & (yar<0),thetaPos+2*np.pi,thetaPos)

        sortedTheta,sortedR=(list(t) for t in zip(*sorted(zip(thetaPos,radialPos))))
        #sortedR=medfilt(sortedR,7)


        fp=np.interp(self.thetaLin,sortedTheta,sortedR,period=2*np.pi)

        return np.array(fp)

    def calibrate(self,image):
        """Function used to get the stationary signal, ie to calibrate"""

        border=self.getBorder(image,1)

        self.calibSig=self.getSignal(border)

        return

    def fourier(self):
        """Calculates the 2D FFT and shifted FFT"""

        if self.signal.shape[0]%2==0:
            self.signal=self.signal[0:-1,:]

        self.numFrames=len(self.signal)
        
        self.fft=np.fft.fft2(self.signal)

        self.fftShifted=np.fft.fftshift(self.fft)

        self.fft=self.fft
        self.fftShifted=np.abs(self.fftShifted)

        self.K2max, self.K1max=self.signal.shape
        self.K2max, self.K1max=int((self.K2max-1)/2),int((self.K1max-1)/2)

        self.kLin=(np.arange(-self.K1max,self.K1max+1))/self.r_phys
        self.kThLin=(np.arange(-self.K1max,self.K1max+1))
        self.omLin=np.arange(-self.K2max,self.K2max+1)*2*np.pi*self.fps/self.numFrames

        self.rngK=range(int(self.D1),int(2*self.K1max-self.D1+1))
        self.rngOm=range(int(self.D2),int(2*self.K2max-self.D2+1))

        self.kLinCrop=self.kLin[self.rngK]
        self.omLinCrop=self.omLin[self.rngOm]
        self.kThLinCrop=self.kThLin[self.rngK]

        self.fftCrop=self.fftShifted[self.rngOm,:][:,self.rngK]

        return

    def drawSpectrum(self):
        """Draws the spectrum cropped around 0,0"""

        figFft,axFft=plt.subplots(figsize=(12,9))

        four=axFft.pcolormesh(self.kLinCrop,self.omLinCrop,np.log(self.fftCrop),vmin=self.cmin,vmax=self.cmax,cmap=parula_map,shading="auto",rasterized=True)
        four.set_edgecolors("face")
        plt.xlabel("k [m$^{-1}$]")
        plt.ylabel("$\omega$ [rad$\cdot$s$^{-1}$]")
        plt.colorbar(four,ax=axFft)
        
        return
    
    def drawTheta(self):
        """Draws the spectrum in dependance on the angular wavenumber"""
        
        plt.subplots(figsize=(12,9))
        
        l1=int(len(self.kThLinCrop+2)/2)
        l2=int(len(self.omLinCrop-1)/2)

        four=plt.pcolormesh(self.kThLinCrop[l1:],self.omLinCrop[l2:],np.log(self.fftCrop[l2:,l1:]),vmin=self.cmin,vmax=self.cmax,cmap=parula_map,shading="auto",rasterized=True)

        four.set_edgecolors("face")

        plt.xlabel(r"k$_\theta$",fontsize=45)
        plt.ylabel("$\omega$ [rad$\cdot$s$^{-1}$]",fontsize=45)
        plt.tick_params(labelsize=30)

        cbar=plt.colorbar(format="%1.1f")
        cbar.set_label(r"log$|\tilde{\eta}(k_\theta,\omega)|$",fontsize=30)
        cbar.ax.tick_params(labelsize=30)
        return


    def dispersion(self,k):
        """Returns the dispersion relation for given k values"""
        
        temp=(self.g_eff*k+self.sigma_eff*k**3/self.rho)*np.tanh(self.c_eff**2*k/self.g_eff)
        
        return np.sqrt(temp)
    
    def dispersionTh(self,k,W,R):
        """Returns the angular dispersion relation for given k values"""
        
        temp=(self.g_eff*k/self.r_phys+self.sigma_eff*k**3/self.rho/self.r_phys**3)*np.tanh(k*W*self.r_phys/R**2/2)

        return np.sqrt(temp)


    def plotDispersion(self,N):
        """Plots the dispersion relation omega(k/N)*N"""

        self.fitDis=self.dispersion(self.kLinCrop/N)*N
        plt.plot(self.kLinCrop[self.kA:-self.kA+1],self.fitDis[self.kA:-self.kA+1], "--r",linewidth=0.7)

        return

    def plotDispersionTh(self,N):
        """Plots the angular dispersion relation omega(k/N)*N"""

        self.fitDisTh=self.dispersionTh(self.kThLinCrop/N)*N
        plt.plot(self.kThLinCrop[self.kA:-self.kA+1],self.fitDisTh[self.kA:-self.kA+1], "--r",linewidth=0.7)

        return


    def plotLinear(self,slope):
        """Plots a linear curve with the given slope"""

        plt.plot(self.kLinCrop,self.kLinCrop*slope,"--r",linewidth=0.7)

        return

    
    
    def saveData(self):
        """ Saves the signal, r_phys, fps, sigma and c into a .npy file with the same name as the video
        in the same location"""

        name=self.videoLoc.split(".mp4")
        name=name[0]
        np.save(name,self.signal)

        return

    def loadData(self):
        """Loads the saved signal etc. assuming that it still has the same name as the video and is located
        in the same folder"""
        
        name=self.videoLoc.split(".mp4")
        name=name[0]

        if(exists(name+".npy")):
            data=np.load(name+".npy",allow_pickle=True)
        else:
            data=np.load(name+"-signal.npy",allow_pickle=True)

        if len(data)>10:
            self.signal=data
        else:
            self.signal=data[0]
        
        return
    

    def filter(self):
        
        self.filtered = np.zeros((len(self.signal),3762))
        
        for i in range(len(self.signal)):
            
            self.filtered[i] = np.convolve(self.signal[i], np.ones(40)/40, mode='valid')
            
        
    def getMin(self):

        self.mins = np.zeros((len(self.filtered),2))
        n2 = np.arange(len(self.filtered[0]))
        
        for i in range(len(self.filtered)):
            
            nxs = argrelextrema(self.filtered[i],np.less,order=5)
            
            arg = n2[nxs][np.argpartition(self.filtered[i][nxs],2)]

            if (len(arg[:2]))==2:
            
                self.mins[i]=arg[:2]

        self.mins = self.mins.astype(int)

        
    def sech(self,x, A, x0, sigma):
        
        return A/np.cosh((x-x0)/sigma)**2

    def getAmp(self):

        self.amps = np.zeros(len(self.filtered))
        self.deltas = np.zeros(len(self.filtered))

        plt.figure()
        
        for i in range(len(self.filtered)):

            x = np.linspace(-1,1,1000)*np.pi*500/len(self.filtered[i])
            
            y = self.filtered[i,self.mins[i,1]-500:self.mins[i,1]+500]

            parameters, covariance = curve_fit(self.sech, x, y)

            h = parameters[0]
            a = parameters[1]
            x0 = parameters[2]
            delta = parameters[3]

            self.amps[i] = a
            self.deltas[i] = delta

            plt.clf()
            plt.plot(x,y)
            plt.plot(x,self.sech(x,h,a,x0,delta))
            
            plt.pause(0.05)
            

    def play(self):
        
        plt.figure()
        
        n2 = np.arange(len(self.filtered[0]))
        
        try: 
            for i in range(len(self.filtered)):
                plt.clf()
                plt.plot(self.filtered[i])
                #plt.ylim([-28,10])
                plt.scatter(n2[self.mins[i][0]],self.filtered[i][self.mins[i][0]])
                plt.scatter(n2[self.mins[i][1]],self.filtered[i][self.mins[i][1]])
                plt.pause(0.05)

        except KeyboardInterrupt:
            pass


    def fit(self,n):

        plt.figure()
        
        plt.plot(sin.signal[n])

        coord = plt.ginput(2)

        dn = (coord[0][0]-coord[1][0])

        n1 = int(coord[0][0])
        n2 = int(coord[1][0])
        
        th = np.linspace(-dn/2,dn/2,n2-n1)*2*np.pi/3801

        dy = 0.5*(coord[0][1]+coord[1][1])

        plt.clf()

        plt.plot(th,sin.signal[n][n1:n2]-dy)

        parameters, covariance = curve_fit(self.sech, th, sin.signal[n][n1:n2]-dy)

        A = parameters[0]
        x0 = parameters[1]
        sigma = parameters[2]
        
        plt.plot(th,self.sech(th, A, x0, sigma))

        print("Delta: ", sigma)
        print("A: ", A)

    def plot(self,n):
        plt.figure()
        
        plt.plot(sin.signal[n])




fname = "/data/torus_wt/09_23/11_09/Basler_acA2040-120um__23597830__20230712_050657186.mp4"


sin = torus(fname,1,550,750)

print("Opening video...")
sin.openVideo()

sin.seuil = 30
sin.showVideo = 0

print("Border detection started")

sin.runVideo()

print()
print("Saving...")
sin.saveData()
print()
print("Done!")
