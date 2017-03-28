#include <time.h>
#include <math.h>
#include "SDL2/SDL.h"
#include "SDL2/SDL_ttf.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include "3D.h"
#include "draw.h"
#include "touch.h"
#include "objloader.h"
#include <windows.h>
#include <inttypes.h>
#include "rs232.h"
#define PORT 19

int fps=0,fpss;
bool done = false;
bool drawdone = false;
pointR3D zero= {0,0,0};
extern char sharebuff[100];

SDL_GameController *controller = NULL;
bool haspad=false;
char buff[1000];

char* getlinecomp(int c)
{
    int count=0;
    char * fp=buff;
    while(*(fp+count-1)!='\n')
    {
        if(RS232_PollComport(c,(unsigned char*)(fp+count),1))
            count++;
    }
    *(fp+count)=0;

    return buff;
}
int16_t buff_i[6];
float buff_f[6];
point3D p3d= {.x=0,.y=0,.z=20};
pointR3D r3d= {.rx=0,.ry=0,.rz=0};
point3D box[8];

void * pth1(void* args)
{
//   char str[30];
//    char str2[30];
    box[0]= {150,0,0};
    box[1]= {0,150,0};
    box[2]= {0,0,-150};
    int objsize;
    ld_list vet;
    ld_list face;
    ld_list tri;
    point3D p1,p2,p3;
    uint8_t * ldbuff=(uint8_t*)malloc(100000000);
    if((objsize=OL_Read("objloader.obj",ldbuff,100000000))==-1)
    {
        return NULL;
    }
    OL_Load(ldbuff,100.0,1,objsize,&vet,&face);
    OL_GetTriangle(&vet,&face,&tri);

    while(!done)
    {
        //usleep(200);
        if(done)
            break;
        clrscreen();
        clrdeepbuff();
        for(uint16_t i=1; i<tri.size; i++)
        {
            OL_Seek(&tri,i);
            p1=setpoint(((node_tri_t*)(tri.p_now))->data.p1.x,((node_tri_t*)(tri.p_now))->data.p1.y,100-((node_tri_t*)(tri.p_now))->data.p1.z,0xffffff);
            p2=setpoint(((node_tri_t*)(tri.p_now))->data.p2.x,((node_tri_t*)(tri.p_now))->data.p2.y,100-((node_tri_t*)(tri.p_now))->data.p2.z,0xffffff);
            p3=setpoint(((node_tri_t*)(tri.p_now))->data.p3.x,((node_tri_t*)(tri.p_now))->data.p3.y,100-((node_tri_t*)(tri.p_now))->data.p3.z,0xffffff);
            //fillthr3D(p1,p2,p3,0xffffff);
            drawline3D(p1,p2,0xffffff,0);
            drawline3D(p3,p2,0xffffff,0);
            drawline3D(p1,p3,0xffffff,0);
            drawline3D(p1,p2,0xffffff,1);
            drawline3D(p3,p2,0xffffff,1);
            drawline3D(p1,p3,0xffffff,1);
        }

        //sprintf(str,"FPS:%d",fpss);
        //sprintf(str2,"Angle of cube:%f",R2D(tet.getangle()));
        //sprintf(sharebuff,"Number of Finger(s):%d",getfingernum());
        //drawstring(0,0,(unsigned char*)str,0xffffff);
        //drawstring(0,10,(unsigned char*)str2,0xffffff);
        drawstring(0,20,(unsigned char*)sharebuff,0xffffff);
        // drawstring(0,30,(unsigned char*)tet.getmsg(),0xffffff);
        //Drawbox(0,0,0,10,r3d,0xfffff);
        //Drawbox(p3d.x,p3d.y,p3d.z,3,zero,0xfffff);
        Drawbox(box[0].x,box[0].y,box[0].z,3,zero,0xfffff);
        for(int i=1; i<3; i++)
            Drawbox(box[i].x,box[i].y,box[i].z,3,zero,0xfff);
        //sprintf(sharebuff,"%3.5f,%3.5f,%3.5f",r3d.rx,r3d.ry,r3d.rz);
        /* for(int i=0;i<4;i++)
         {
         sprintf(str2,"Key%d:%d",i,axis[i]);

         drawstring(0,40+i*10,(unsigned char*)str2,0xffffff);

         }*/

        updatescreen();
        fps++;
    }
    drawdone=true;
    return NULL;
}

void * pth2(void* args)///-----------------------------------这是fps计数线程
{
    while(1)
    {

        fps=0;
        for(int i=0; i<10; i++)
            usleep(100000);
        fpss=fps;
    }
    return NULL;
}
bool doclear=false;
void * pth3(void* args)
{
    while(1)
    {
        usleep(10000);
        uint32_t key=getkey();
        if(key&0x01)
            doclear=true;

        //movecam(0,0,1);
        /*      if(key&0x02)
                  movecam(0,0,-1);
              if(key&0x04)
                  movecam(1,0,0);
              if(key&0x08)
                  movecam(-1,0,0);
              if(key&0x10)
                  movecam(0,-1,0);
              if(key&0x20)
                  movecam(0,1,0);*/
    }

    return NULL;
}




pointR3D R3D_Axis(pointR3D& input,pointR3D& axis,float angle)
{
    pointR3D temp;
    /* temp.rx=(axis.rx*axis.rx+(1-axis.rx*axis.rx)*cos_)*input.rx+\
             (axis.rx*axis.ry*(1-cos_)+axis.rz*sinf(angle))*input.ry+\
             (axis.rx*axis.rz*(1-cos_)-axis.ry*sinf(angle))*input.rz;
     temp.ry=(axis.rx*axis.ry*(1-cos_)-axis.rz*sinf(angle))*input.rx+\
             (axis.ry*axis.ry+(1-axis.ry*axis.ry)*cos_)*input.ry+\
             (axis.rz*axis.ry*(1-cos_)+axis.rx*sinf(angle))*input.rz;
     temp.rz=(axis.rx*axis.rz*(1-cos_)+axis.ry*sinf(angle))*input.rx+\
             (axis.rz*axis.ry*(1-cos_)-axis.rx*sinf(angle))*input.ry+\
             (axis.rz*axis.rz+(1-axis.rz*axis.rz)*cos_)*input.rz;*/
    float sin_=sinf(angle),cos_=cosf(angle);
    temp.rx=(axis.rx*axis.rx+(1-axis.rx*axis.rx)*cos_)*input.rx+\
            (axis.rx*axis.ry*(1-cos_)-axis.rz*sin_)*input.ry+\
            (axis.rx*axis.rz*(1-cos_)+axis.ry*sin_)*input.rz;
    temp.ry=(axis.rx*axis.ry*(1-cos_)+axis.rz*sin_)*input.rx+\
            (axis.ry*axis.ry+(1-axis.ry*axis.ry)*cos_)*input.ry+\
            (axis.rz*axis.ry*(1-cos_)-axis.rx*sin_)*input.rz;
    temp.rz=(axis.rx*axis.rz*(1-cos_)-axis.ry*sin_)*input.rx+\
            (axis.rz*axis.ry*(1-cos_)+axis.rx*sin_)*input.ry+\
            (axis.rz*axis.rz+(1-axis.rz*axis.rz)*cos_)*input.rz;
    return temp;
}

void * pth4(void* args)
{
    pointR3D axis[3];
    pointR3D axis_buff[3];
    axis[0]= {1.0,0,0};
    axis[1]= {0,1.0,0};
    axis[2]= {0,0,1.0};
    pointR3D cali_g= {0.000050,0.000009,-0.000021};
    /*point3D Rest= {0,0,0};
    auto get = [ cali_g]()
    {
        if(sscanf(getlinecomp(PORT),"%hd/%hd/%hd,%hd/%hd/%hd",&buff_i[0],&buff_i[1],&buff_i[2],&buff_i[3],&buff_i[4],&buff_i[5])==6)
        {
            buff_f[0]=buff_i[0]/32768.0f*16.0;
            buff_f[1]=buff_i[1]/32768.0f*16.0;
            buff_f[2]=buff_i[2]/32768.0f*16.0;
            buff_f[3]=D2R(buff_i[3]/32768.0f)-cali_g.rx;///drx
            buff_f[5]=-D2R(buff_i[4]/32768.0f)-cali_g.ry;///dry
            buff_f[4]=D2R(buff_i[5]/32768.0f)-cali_g.rz;///drz
        }
    };

    for(int i=0; i<100; i++)
    {
        get();
        Rest.x-=buff_f[1];
        Rest.y+=buff_f[2];
        Rest.z+=buff_f[0];
        cali_g.rz+=buff_f[3];
        cali_g.rx+=buff_f[4];
        cali_g.ry+=buff_f[5];
    }
    Rest.x/=100.0;
    Rest.y/=100.0;
    Rest.z/=100.0;
    box[0].x=Rest.x*15;
    box[0].y=Rest.y*15;
    box[0].z=Rest.z*15;
    cali_g.rx/=100.0;
    cali_g.rx+=0.000070;
    cali_g.ry/=100.0;
    cali_g.ry-=0.0000001;
    cali_g.rz/=100.0;
    cali_g.rz-=0.000070;
    sprintf(sharebuff,"%f,%f,%f",cali_g.rx,cali_g.ry,cali_g.rz);
    box[4].x=Rest.x;
    box[4].y=Rest.y;
    box[4].z=Rest.z;*/
    while(!done)
    {
        //usleep(2);
        static    auto get = [ cali_g]()
        {
            if(sscanf(getlinecomp(PORT),"%hd/%hd/%hd,%hd/%hd/%hd",&buff_i[0],&buff_i[1],&buff_i[2],&buff_i[3],&buff_i[4],&buff_i[5])==6)
            {
                buff_f[0]=buff_i[0]/32768.0f*16.0;
                buff_f[1]=buff_i[1]/32768.0f*16.0;
                buff_f[2]=buff_i[2]/32768.0f*16.0;
                buff_f[3]=D2R(buff_i[3]/32768.0f)-cali_g.rx;///drx
                buff_f[5]=-D2R(buff_i[4]/32768.0f)-cali_g.ry;///dry
                buff_f[4]=D2R(buff_i[5]/32768.0f)-cali_g.rz;///drz
            }
        };

        get();
        pointR3D Rgyro= {buff_f[4],-buff_f[5],-buff_f[3]};
        for(int i=0; i<3; i++)
        {
            axis_buff[i]=axis[i];
        }
        axis_buff[1]=R3D_Axis(axis_buff[1],axis[0],Rgyro.rx);
        axis_buff[2]=R3D_Axis(axis_buff[2],axis[0],Rgyro.rx);
        axis_buff[0]=R3D_Axis(axis_buff[0],axis[1],Rgyro.ry);
        axis_buff[2]=R3D_Axis(axis_buff[2],axis[1],Rgyro.ry);
        axis_buff[0]=R3D_Axis(axis_buff[0],axis[2],Rgyro.rz);
        axis_buff[1]=R3D_Axis(axis_buff[1],axis[2],Rgyro.rz);
        for(int i=0; i<3; i++)
        {
            axis[i]=axis_buff[i];
        }
        if(doclear)
        {
            axis[0]= {1.0,0,0};
            axis[1]= {0,1.0,0};
            axis[2]= {0,0,1.0};
            doclear=false;
        }
        setvraxis(axis[0],axis[1],axis[2]);
    }
    return NULL;
}

void * pth5(void* args)
{
    while(!done)
    {
        usleep(5000);
    }
    return NULL;
}

int main(int argc, char *argv[])
{
    if(RS232_OpenComport(PORT,115200))
        sprintf(sharebuff,"Open COM Port Failed!");
    pthread_t thread,thread2,thread3,thread4;
    setdisplay(1280,720);
    initSDL();
    for (int i = 0; i < SDL_NumJoysticks(); ++i)
        if (SDL_IsGameController(i))
        {
            controller = SDL_GameControllerOpen(i);
            if (controller)
            {
                haspad=true;
                break;

            }
        }
    setcam(setpoint(0,0,0,0),setpointR(0,0,0));
    seteye(setpoint(0,0,-1000,0));
    pthread_create(&thread,NULL,pth1,NULL);
    pthread_create(&thread2,NULL,pth2,NULL);
    pthread_create(&thread3,NULL,pth3,NULL);
    pthread_create(&thread4,NULL,pth4,NULL);
    SDLloop(&drawdone,&done);
    exit(0);
}






