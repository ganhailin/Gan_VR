#include "SDL2/SDL.h"
#include <math.h>
#include <stdlib.h>
#include "3D.h"
#include "draw.h"
#include "ASICC.h"
#include <string.h>
#include <stdio.h>
#define _DEEPBUFF

#define TEST_
#define xmax 1920
#define ymax 1080
static float deepbuff[xmax][ymax];

pointR3D rd0,rd1;
pointR3D camR0;
point3D cam0;
point3D eye0;
pointR3D xaxis;
pointR3D yaxis;
pointR3D zaxis;


void Swap(int &x, int &y)
{
    int temp = y;
    y = x;
    x = temp;
}

inline void Swapf(float &x, float &y)
{
    float temp = y;
    y = x;
    x = temp;
}

inline float vet_dot(pointR3D a,pointR3D b)
{
    return a.rx*b.rx+a.ry*b.ry+a.rz*b.rz;
}
inline float vet_dot(point3D a,pointR3D b)
{
    return a.x*b.rx+a.y*b.ry+a.z*b.rz;
}
inline float abs(pointR3D a)
{
    return sqrtf(a.rx*a.rx+a.ry*a.ry+a.rz*a.rz);
}
pointR3D la={0,5.f/180.f*3.1415926535897,0},ra={0,-5.f/180.f*3.1415926535897,0};

int point3to3(point3D p3d, point3D * p2d,int s)
{
    point3D temp,temp2;
    temp2.x = -p3d.x - cam0.x;
    temp2.y = p3d.y - cam0.y;
    temp2.z = p3d.z - cam0.z;
    temp.x=-vet_dot(temp2,xaxis);
    temp.y=vet_dot(temp2,yaxis);
    temp.z=vet_dot(temp2,zaxis);
    temp.x-=s?-10:10;

   // if (!(camR0.rx == 0 && camR0.ry == 0 && camR0.rz == 0))
   if(s)
        R3D(&temp, &la);
        else
            R3D(&temp, &ra);
    p2d->x =  (temp.x - eye0.x) * ((eye0.z) / temp.z);
    p2d->y =  (temp.y - eye0.y) * ((eye0.z) / temp.z);
    p2d->z = temp.z;
    p2d->c = p3d.c;
    p2d->x+=(float)xrad / 4;
    p2d->y+=(float)yrad / 2;
    return temp.z > 1 ? 1 : 0;
}

void setvraxis(pointR3D &x,pointR3D &y,pointR3D &z)
{
    xaxis=x;
    yaxis=y;
    zaxis=z;
}
/*
void DrawPixel3D(int x, int y, float z, int c)
{
    if (x >= 0 && x < xrad && y >= 0 && y < yrad)
        if (z <= deepbuff[x][y])
        {
            deepbuff[x][y] = z;
            Draw_Pixel(x, y, c);
        }
}*/
void DrawPixel3D(int x, int y, float z, int c,int s)
{
    if(s)
    {
        if (x >= 0&& x < xrad/2 && y >= 0 && y < yrad)
            if (z <= deepbuff[x+xrad/2][y])
            {
                deepbuff[x+xrad/2][y] = z;
                Draw_Pixel(x+xrad/2, y, c);
            }

    }
    else
    {
        if (x >= 0 && x < xrad/2 && y >= 0 && y < yrad)
            if (z <= deepbuff[x][y])
            {
                deepbuff[x][y] = z;
                Draw_Pixel(x, y, c);
            }

    }
}

void clrdeepbuff(void)
{
    for (int x = 0; x < xrad; x++)
        for (int y = 0; y < yrad; y++)
            deepbuff[x][y] = 4.2E38;
}

void Draw_Simple_Line(float x1, float x2, float z1, float z2, float y, int color,int s)
{
    Drawline3D(x1, y, z1, x2, y, z2, color,s);
}

int Draw_Top_Trangle(float x0, float y0, float z0,
                     float x1, float y1, float z1,
                     float x2, float y2, float z2, int color,int s)
{

    if (y0 == y1)
    {
    }
    else if (y0 == y2)
    {
        Swapf(x2, x1);
        Swapf(y2, y1);
        Swapf(z2, z1);
    }
    else if (y1 == y2)
    {
        Swapf(x0, x2);
        Swapf(y0, y2);
        Swapf(z0, z2);
    }
    else
    {
        return 1;
    }

    if (x1 < x0)
    {
        Swapf(x1, x0);
        Swapf(y1, y0);
        Swapf(z1, z0);
    }
    else if (x1 == x0)
    {
        return 1;
    }

    float dxy_left = (x2 - x0)/(y2 - y0);
    float dxy_right = (x1 - x2) / (y1 - y2);
    float dzy_left = (z2 - z0)/ (y2 - y0);
    float dzy_right = (z2 - z1)/ (y2 - y1);

    float xs = x0, xe = x1, zs = z0, ze = z1;
    for (float y = y0; y <= y2; y++)
    {
        Draw_Simple_Line(x0+(y-y0)*dxy_left, x1+(y-y0)*dxy_right, z0+(y-y0)*dzy_left,z1+(y-y0)*dzy_right, y, color,s);
        zs += dzy_left;
        ze += dzy_right;

        xs += dxy_left;
        xe += dxy_right;
    }
    return 0;
}
int Draw_Bottom_Trangle(float x0, float y0, float z0,
                        float x1, float y1, float z1,
                        float x2, float y2, float z2, int color,int s)
{
    if (y2 == y1)
    {
    }
    else if (y2 == y0)
    {
        Swapf(x0, x1);
        Swapf(y0, y1);
        Swapf(z0, z1);
    }
    else if (y0 == y1)
    {
        Swapf(x0, x2);
        Swapf(y0, y2);
        Swapf(z0, z2);
    }
    else
    {
        return 1;
    }

    if (x1 < x2)
    {
        Swapf(x1, x2);
        Swapf(y1, y2);
        Swapf(z1, z2);
    }
    else if (x1 == x2)
    {
        return 1;
    }

    float dxy_left = (x2 - x0)/(y2 - y0);
    float dxy_right = (x1 - x0)/(y1 - y0);
    float dzy_left = (z2 - z0)/(y2 - y0);
    float dzy_right = (z1 - z0)/(y1 - y0);

    float xs = x0, xe = x0, zs = z0, ze = z0;
    for (float y = y0; y <= y2; y++)
    {
        Draw_Simple_Line(xs, xe, zs,ze, y, color,s);
        zs += dzy_left;
        ze += dzy_right;
        xs += dxy_left;
        xe += dxy_right;
    }
    return 0;
}

int Draw_Trangle_2D(float x0, float y0, float z0,
                    float x1, float y1, float z1,
                    float x2, float y2, float z2,
                    int color,int s)
{
    if (!((x0 >= -100 && x0 < xrad+100 && y0 >= -100 && y0 < yrad+100)&&(x1 >= -100 && x1 < xrad+100 && y1 >= -100 && y1 < yrad+100)&&(x2 >= -100 && x2 < xrad+100 && y2 >= -100 && y2 < yrad+100)))return 0;

    if ((x0 == x1 && x1 == x2) || (y0 == y1 && y1 == y2))
    {
        return 1;
    }

    if (y0 > y1)
    {
        Swapf(x0, x1);
        Swapf(y0, y1);
        Swapf(z0, z1);
    }

    if (y0 > y2)
    {
        Swapf(x0, x2);
        Swapf(y0, y2);
        Swapf(z0, z2);
    }

    if (y1 > y2)
    {
        Swapf(y1, y2);
        Swapf(x1, x2);
        Swapf(z1, z2);
    }

    if (y0 == y1)
    {
        Draw_Top_Trangle(x0, y0, z0, x1, y1, z1, x2, y2, z2, color,s);
    }
    else if (y1 == y2)
    {
        Draw_Bottom_Trangle(x0, y0, z0, x1, y1, z1, x2, y2, z2, color,s);
    }
    else
    {
        float new_x = x0 + (y1 - y0) * (x2 - x0) / (y2 - y0);
        float new_z =z0+(y1 - y0)*(z2-z0)  / (y2 - y0);
        Draw_Bottom_Trangle(x0, y0, z0, new_x, y1, new_z, x1, y1, z1, color,s);
        Draw_Top_Trangle(new_x, y1, new_z, x1, y1, z1, x2, y2, z2, color,s);
    }


    return 0;
}



int Drawline3D(int x1, int y1, float z1, int x2, int y2, float z2, int color,int s)
{
    if (!((x1 >= -300 && x1 < xrad+300 && y1 >= -300 && y1 < yrad+300)&&(x2 >= -300 && x2 < xrad+300 && y2 >= -300 && y2 < yrad+300)))return 0;
    int dx = x2 - x1, dy = y2 - y1;
    float dz = z2 - z1;
    float dmax = abs(dx) > abs(dy) ? abs(dx) : abs(dy);
    for (int i = 0; i < dmax; i++)
    {
        DrawPixel3D(x1 + dx * i / dmax, y1 + dy * i / dmax, z1 + dz *1.0*( i / dmax),
                    color,s);

    }
    DrawPixel3D(x2, y2, z2, color,s);
    return 0;
}
void fillthr3D(point3D p1, point3D p2, point3D p3, int color)
{
    point3D b1=p1,b2=p2,b3=p3;
    if (point3to3(b1, &p1,1) && point3to3(b2, &p2,1) && point3to3(b3, &p3,1))
    {
        Draw_Trangle_2D(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,p3.x,p3.y,p3.z,color,1);
    }
    if (point3to3(b1, &p1,0) && point3to3(b2, &p2,0) && point3to3(b3, &p3,0))
    {
        Draw_Trangle_2D(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,p3.x,p3.y,p3.z,color,0);
    }
}



void DrawPixel(int x, int y, int c)
{
    if (x >= 0 && x < xrad && y >= 0 && y < yrad)
        Draw_Pixel(x, y, c);
}


int Drawline(int x1, int y1, int x2, int y2, int color)
{
    int y, x;

    if (y1 == y2)
    {
        y = y1;
        if (x1 >= x2)
            for (x = x2; x < x1; x++)
                DrawPixel(x, y, color);
        else
            for (x = x1; x < x2; x++)
                DrawPixel(x, y, color);
    }
    if (x1 == x2)
    {
        x = x1;
        if (y1 <= y2)
            for (y = y1; y < y2; y++)
                DrawPixel(x, y, color);
        else
            for (y = y2; y < y1; y++)
                DrawPixel(x, y, color);


    }
    else
    {
        int dx = abs(x1 - x2), dy = abs(y1 - y2);

        if (dy > dx)
        {
            if (y1 < y2)
            {
                for (y = y1; y <= y2; y++)
                    DrawPixel(x1 + (x2 - x1) * (y - y1) / dy, y, color);
            }
            else
                Drawline(x2, y2, x1, y1, color);
        }
        else
        {
            if (x1 < x2)
            {
                for (x = x1; x <= x2; x++)
                    DrawPixel(x, y1 + (y2 - y1) * (x - x1) / dx, color);
            }
            else
                Drawline(x2, y2, x1, y1, color);

        }
    }
    return 0;
}


void fillthr(Point p1, Point p2, Point p3, int color)
{

    float dec1, dec2, d1, d2;
    int y;
    Point p4;

    if (p1.y == p2.y)
    {
        dec1 = (float)(((float)(p1.y - p3.y)) / (float)((p1.x - p3.x)));
        d1 = p1.y - dec1 * p1.x;
        dec2 = (float)(((float)(p2.y - p3.y)) / (float)((p2.x - p3.x)));
        d2 = p2.y - dec2 * p2.x;

        if (p3.y > p2.y)
            if (p1.x == p3.x)
                for (y = p2.y; y <= p3.y; y++)
                    Drawline(p3.x, y, (int)((y - d2) / dec2), y, color);
            else if (p2.x == p3.x)
                for (y = p2.y; y <= p3.y; y++)
                    Drawline((int)((y - d1) / dec1), y, p3.x, y, color);
            else
                for (y = p2.y; y <= p3.y; y++)
                    Drawline((int)((y - d1) / dec1), y, (int)((y - d2) / dec2),
                             y, color);
        else if (p1.x == p3.x)
            for (y = p3.y; y <= p2.y; y++)
                Drawline(p3.x, y, (int)((y - d2) / dec2), y, color);
        else if (p2.x == p3.x)
            for (y = p3.y; y <= p2.y; y++)
                Drawline((int)((y - d1) / dec1), y, p3.x, y, color);
        else
            for (y = p3.y; y <= p2.y; y++)
                Drawline((y - d1) / dec1, y, (y - d2) / dec2, y, color);
    }
    else
    {
        if (p2.y == p3.y)
            fillthr(p2, p3, p1, color);
        else
        {
            if (p1.y == p3.y)
                fillthr(p1, p3, p2, color);
            else
            {
                if (((p2.y < p1.y) && (p2.y > p3.y))|| ((p2.y > p1.y) && (p2.y < p3.y)))
                {
                    p4.y = p2.y;
                    dec1 =
                        (float)(((float)(p1.y - p3.y)) /
                                ((float)(p1.x - p3.x)));
                    d1 = (float)(((float)p1.y) - ((float)(dec1 * p1.x)));
                    p4.x = (p4.y - d1) / dec1;
                    fillthr(p1, p4, p2, color);
                    fillthr(p3, p4, p2, color);
                }
                else
                {
                    if (((p1.y < p2.y) && (p1.y > p3.y))|| ((p2.y < p1.y) && (p1.y < p3.y)))
                        fillthr(p3, p1, p2, color);
                    else
                    {
                        if (((p3.y < p1.y) && (p3.y > p2.y)) || ((p3.y > p1.y) && (p3.y < p2.y)))
                            fillthr(p1, p3, p2, color);
                    }
                }
            }
        }
    }
}


void R3Dx(point3D * p3d, float rx)
{
    double y, z;
    y = cos(rx) * p3d->y - sin(rx) * p3d->z;
    z = cos(rx) * p3d->z + sin(rx) * p3d->y;
    p3d->y = y;
    p3d->z = z;
}

void R3Dy(point3D * p3d, float ry)
{
    double x, z;
    x = cos(ry) * p3d->x + sin(ry) * p3d->z;
    z = -sin(ry) * p3d->x + cos(ry) * p3d->z;
    p3d->x = x;
    p3d->z = z;
}

void R3Dz(point3D * p3d, float rz)
{
    double x, y;
    x = cos(rz) * p3d->x - sin(rz) * p3d->y;
    y = cos(rz) * p3d->y + sin(rz) * p3d->x;
    p3d->x = x;
    p3d->y = y;
}

void R3D(point3D * p3d, pointR3D * p3dr)
{

    float cosrx = cosf(p3dr->rx),
          cosry = cosf(-p3dr->ry),
          cosrz = cosf(-p3dr->rz),
          sinrx = sinf(p3dr->rx),
          sinry = sinf(-p3dr->ry),
          sinrz = sinf(-p3dr->rz);
    float x =
        (cosry * cosrz) * p3d->x +
        (-cosry * sinrz) * p3d->y + (sinry) * p3d->z;
    float y =
        (cosrz * sinrx * sinry + cosrx * sinrz) * p3d->x +
        (cosrx * cosrz - sinrx * sinry * sinrz) * p3d->y +
        (-cosry * sinrx) * p3d->z;
    float z =
        (-cosrx * cosrz * sinry + sinrx * sinrz) * p3d->x +
        (cosrz * sinrx + cosrx * sinry * sinrz) * p3d->y +
        (cosrx * cosry) * p3d->z;
    p3d->x = x;
    p3d->y = y;
    p3d->z = z;
}

void R3DM(point3D * p3d, M_R3D rm, int count)
{
    for (int i = 0; i < count; i++)
    {
        float x = rm.x1 * p3d[i].x + p3d[i].y * rm.x2 + rm.x3 * p3d[i].z;
        float y = rm.y1 * p3d[i].x + p3d[i].y * rm.y2 + rm.y3 * p3d[i].z;
        float z = rm.z1 * p3d[i].x + p3d[i].y * rm.z2 + rm.z3 * p3d[i].z;
        p3d[i].x = x;
        p3d[i].y = y;
        p3d[i].z = z;
    }
}


void getrm(M_R3D * rm, pointR3D * p3dr)
{

    float cosrx = cosf(p3dr->rx),
          cosry = cosf(-p3dr->ry),
          cosrz = cosf(-p3dr->rz),
          sinrx = sinf(p3dr->rx),
          sinry = sinf(-p3dr->ry),
          sinrz = sinf(-p3dr->rz);

    rm->x1 = (cosry * cosrz);
    rm->x2 = (-cosry * sinrz);
    rm->x3 = (sinry);
    rm->y1 = (cosrz * sinrx * sinry + cosrx * sinrz);
    rm->y2 = (cosrx * cosrz - sinrx * sinry * sinrz);
    rm->y3 = (-cosry * sinrx);
    rm->z1 = (-cosrx * cosrz * sinry + sinrx * sinrz);
    rm->z2 = (cosrz * sinrx + cosrx * sinry * sinrz);
    rm->z3 = (cosrx * cosry);
}

int point3to2(point3D p3d, point * p2d)
{
    point3D temp;
    temp.x = -p3d.x - cam0.x;
    temp.y = p3d.y - cam0.y;
    temp.z = p3d.z - cam0.z;
    if (!(camR0.rx == 0 && camR0.ry == 0 && camR0.rz == 0))
        R3D(&temp, &camR0);
    p2d->x = xrad / 2 + (temp.x - eye0.x) * ((eye0.z) / temp.z);
    p2d->y = yrad / 2 + (temp.y - eye0.y) * ((eye0.z) / temp.z);
    p2d->c = p3d.c;
    return temp.z > 1 ? 1 : 0;
}


void drawline3D(point3D p1, point3D p2, int color,int s)
{
    if (point3to3(p1, &p1,s))
        if (point3to3(p2, &p2,s))
        {
            Drawline3D(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, color,s);
        }
}

void Drawbox(float x, float y, float z, float r, pointR3D rd, int color)
{
    point3D box[8];
    M_R3D rm;
    box[0].x = x-r;
    box[0].y = y-r;
    box[0].z = z-r;
    box[1].x = x+r;
    box[1].y = y-r;
    box[1].z = z-r;
    box[2].x = x-r;
    box[2].y = y+r;
    box[2].z = z-r;
    box[3].x = x+r;
    box[3].y = y+r;
    box[3].z = z-r;
    box[4].x = x-r;
    box[4].y = y-r;
    box[4].z = z+r;
    box[5].x = x+r;
    box[5].y = y-r;
    box[5].z = z+r;
    box[6].x = x-r;
    box[6].y = y+r;
    box[6].z = z+r;
    box[7].x = x+r;
    box[7].y = y+r;
    box[7].z = z+r;
    getrm(&rm, &rd);
    R3DM(box, rm, 8);
#ifdef TEST_
//    int color2=(((color>>4)&0xff)*2/3)<<4|(((color>>2)&0xff)*2/3)<<2|((color&0xff)*2/3);
    fillthr3D(box[0],box[1],box[3],color);
    fillthr3D(box[0],box[2],box[3],color);
    fillthr3D(box[0],box[1],box[5],color);
    fillthr3D(box[0],box[4],box[5],color);
    fillthr3D(box[2],box[3],box[7],color);
    fillthr3D(box[2],box[6],box[7],color);
    fillthr3D(box[6],box[7],box[5],color);
    fillthr3D(box[6],box[4],box[5],color);
    fillthr3D(box[0],box[2],box[6],color);
    fillthr3D(box[0],box[4],box[6],color);
    fillthr3D(box[1],box[5],box[7],color);
    fillthr3D(box[1],box[3],box[7],color);/*
    drawline3D(box[4], box[6], color2);
    drawline3D(box[4], box[5], color2);
    drawline3D(box[5], box[7], color2);
    drawline3D(box[6], box[7], color2);
    drawline3D(box[0], box[1], color2);
    drawline3D(box[0], box[2], color2);
    drawline3D(box[0], box[4], color2);
    drawline3D(box[1], box[5], color2);
    drawline3D(box[1], box[3], color2);
    drawline3D(box[2], box[3], color2);
    drawline3D(box[2], box[6], color2);
    drawline3D(box[3], box[7], color2);*/

#else
    drawline3D(box[4], box[6], color);
    drawline3D(box[4], box[5], color);
    drawline3D(box[5], box[7], color);
    drawline3D(box[6], box[7], color);
    drawline3D(box[0], box[1], color);
    drawline3D(box[0], box[2], color);
    drawline3D(box[0], box[4], color);
    drawline3D(box[1], box[5], color);
    drawline3D(box[1], box[3], color);
    drawline3D(box[2], box[3], color);
    drawline3D(box[2], box[6], color);
    drawline3D(box[3], box[7], color);
#endif // TEST_
}
void Drawbox2(float x, float y, float z, float r, pointR3D rd,pointR3D rd0, int color)
{
    point3D box[8];
    M_R3D rm;
    box[0].x = -r;
    box[0].y = -r;
    box[0].z = -r;
    box[1].x = +r;
    box[1].y = -r;
    box[1].z = -r;
    box[2].x = -r;
    box[2].y = +r;
    box[2].z = -r;
    box[3].x = +r;
    box[3].y = +r;
    box[3].z = -r;
    box[4].x = -r;
    box[4].y = -r;
    box[4].z = +r;
    box[5].x = +r;
    box[5].y = -r;
    box[5].z = +r;
    box[6].x = -r;
    box[6].y = +r;
    box[6].z = +r;
    box[7].x = +r;
    box[7].y = +r;
    box[7].z = +r;
    getrm(&rm, &rd);
    R3DM(box, rm, 8);
    box[0].x += x;
    box[0].y += y;
    box[0].z += z;
    box[1].x += x;
    box[1].y += y;
    box[1].z += z;
    box[2].x += x;
    box[2].y += y;
    box[2].z += z;
    box[3].x += x;
    box[3].y += y;
    box[3].z += z;
    box[4].x += x;
    box[4].y += y;
    box[4].z += z;
    box[5].x += x;
    box[5].y += y;
    box[5].z += z;
    box[6].x += x;
    box[6].y += y;
    box[6].z += z;
    box[7].x += x;
    box[7].y += y;
    box[7].z += z;
    getrm(&rm, &rd0);
    R3DM(box, rm, 8);
#ifdef TEST_
    fillthr3D(box[0],box[1],box[3],color);
    fillthr3D(box[0],box[2],box[3],color);
    fillthr3D(box[0],box[1],box[5],color);
    fillthr3D(box[0],box[4],box[5],color);
    fillthr3D(box[2],box[3],box[7],color);
    fillthr3D(box[2],box[6],box[7],color);
    fillthr3D(box[6],box[7],box[5],color);
    fillthr3D(box[6],box[4],box[5],color);
    fillthr3D(box[0],box[2],box[6],color);
    fillthr3D(box[0],box[4],box[6],color);
    fillthr3D(box[1],box[5],box[7],color);
    fillthr3D(box[1],box[3],box[7],color);/*
    int color2=(((color/65535)&0xff)*2/3)<<4|(((color=256)&0xff)*2/3)<<2|((color&0xff)*2/3);
    drawline3D(box[4], box[6], color2);
    drawline3D(box[4], box[5], color2);
    drawline3D(box[5], box[7], color2);
    drawline3D(box[6], box[7], color2);
    drawline3D(box[0], box[1], color2);
    drawline3D(box[0], box[2], color2);
    drawline3D(box[0], box[4], color2);
    drawline3D(box[1], box[5], color2);
    drawline3D(box[1], box[3], color2);
    drawline3D(box[2], box[3], color2);
    drawline3D(box[2], box[6], color2);
    drawline3D(box[3], box[7], color2);*/

#else
    drawline3D(box[4], box[6], color);
    drawline3D(box[4], box[5], color);
    drawline3D(box[5], box[7], color);
    drawline3D(box[6], box[7], color);
    drawline3D(box[0], box[1], color);
    drawline3D(box[0], box[2], color);
    drawline3D(box[0], box[4], color);
    drawline3D(box[1], box[5], color);
    drawline3D(box[1], box[3], color);
    drawline3D(box[2], box[3], color);
    drawline3D(box[2], box[6], color);
    drawline3D(box[3], box[7], color);
#endif // TEST_
}

void DrawboxM(float x, float y, float z, float r, pointR3D rd, int color, int is)
{

    point3D box[8];
    M_R3D rm;
    box[0].x = x-r;
    box[0].y = y-r;
    box[0].z = z-r;
    box[1].x = x+r;
    box[1].y = y-r;
    box[1].z = z-r;
    box[2].x = x-r;
    box[2].y = y+r;
    box[2].z = z-r;
    box[3].x = x+r;
    box[3].y = y+r;
    box[3].z = z-r;
    box[4].x = x-r;
    box[4].y = y-r;
    box[4].z = z+r;
    box[5].x = x+r;
    box[5].y = y-r;
    box[5].z = z+r;
    box[6].x = x-r;
    box[6].y = y+r;
    box[6].z = z+r;
    box[7].x = x+r;
    box[7].y = y+r;
    box[7].z = z+r;
    if (is > 0)
        getrm(&rm, &rd);
    R3DM(box, rm, 8);
#ifdef TEST_

    fillthr3D(box[0],box[1],box[3],color);
    fillthr3D(box[0],box[2],box[3],color);
    fillthr3D(box[0],box[1],box[5],color);
    fillthr3D(box[0],box[4],box[5],color);
    fillthr3D(box[2],box[3],box[7],color);
    fillthr3D(box[2],box[6],box[7],color);
    fillthr3D(box[6],box[7],box[5],color);
    fillthr3D(box[6],box[4],box[5],color);
    fillthr3D(box[0],box[2],box[6],color);
    fillthr3D(box[0],box[4],box[6],color);
    fillthr3D(box[1],box[5],box[7],color);
    fillthr3D(box[1],box[3],box[7],color);
    fillthr3D(box[0],box[1],box[2],color);
    fillthr3D(box[1],box[2],box[3],color);
    fillthr3D(box[0],box[1],box[4],color);
    fillthr3D(box[1],box[4],box[5],color);
    fillthr3D(box[2],box[3],box[6],color);
    fillthr3D(box[3],box[6],box[7],color);
    fillthr3D(box[6],box[7],box[4],color);
    fillthr3D(box[7],box[4],box[5],color);
    fillthr3D(box[0],box[2],box[4],color);
    fillthr3D(box[2],box[4],box[6],color);
    fillthr3D(box[1],box[5],box[3],color);
    fillthr3D(box[5],box[3],box[7],color);/*
    int color2=0;
    drawline3D(box[4], box[6], color2);
    drawline3D(box[4], box[5], color2);
    drawline3D(box[5], box[7], color2);
    drawline3D(box[6], box[7], color2);
    drawline3D(box[0], box[1], color2);
    drawline3D(box[0], box[2], color2);
    drawline3D(box[0], box[4], color2);
    drawline3D(box[1], box[5], color2);
    drawline3D(box[1], box[3], color2);
    drawline3D(box[2], box[3], color2);
    drawline3D(box[2], box[6], color2);
    drawline3D(box[3], box[7], color2);*/

#else
    drawline3D(box[4], box[6], color);
    drawline3D(box[4], box[5], color);
    drawline3D(box[5], box[7], color);
    drawline3D(box[6], box[7], color);
    drawline3D(box[0], box[1], color);
    drawline3D(box[0], box[2], color);
    drawline3D(box[0], box[4], color);
    drawline3D(box[1], box[5], color);
    drawline3D(box[1], box[3], color);
    drawline3D(box[2], box[3], color);
    drawline3D(box[2], box[6], color);
    drawline3D(box[3], box[7], color);
#endif // TEST_

}

void DrawboxM2(float x, float y, float z, float r, pointR3D rd,pointR3D rd2, int color, int is)
{

    point3D box[8];
    M_R3D rm,rm2;
    box[0].x = -r;
    box[0].y = -r;
    box[0].z = -r;
    box[1].x = +r;
    box[1].y = -r;
    box[1].z = -r;
    box[2].x = -r;
    box[2].y = +r;
    box[2].z = -r;
    box[3].x = +r;
    box[3].y = +r;
    box[3].z = -r;
    box[4].x = -r;
    box[4].y = -r;
    box[4].z = +r;
    box[5].x = +r;
    box[5].y = -r;
    box[5].z = +r;
    box[6].x = -r;
    box[6].y = +r;
    box[6].z = +r;
    box[7].x = +r;
    box[7].y = +r;
    box[7].z = +r;
    if (is > 0)
        getrm(&rm, &rd);
    R3DM(box, rm, 8);
    box[0].x += x;
    box[0].y += y;
    box[0].z += z;
    box[1].x += x;
    box[1].y += y;
    box[1].z += z;
    box[2].x += x;
    box[2].y += y;
    box[2].z += z;
    box[3].x += x;
    box[3].y += y;
    box[3].z += z;
    box[4].x += x;
    box[4].y += y;
    box[4].z += z;
    box[5].x += x;
    box[5].y += y;
    box[5].z += z;
    box[6].x += x;
    box[6].y += y;
    box[6].z += z;
    box[7].x += x;
    box[7].y += y;
    box[7].z += z;
    if (is > 0)
        getrm(&rm2, &rd2);
    R3DM(box, rm2, 8);
#ifdef TEST_

    fillthr3D(box[0],box[1],box[3],color);
    fillthr3D(box[0],box[2],box[3],color);
    fillthr3D(box[0],box[1],box[5],color);
    fillthr3D(box[0],box[4],box[5],color);
    fillthr3D(box[2],box[3],box[7],color);
    fillthr3D(box[2],box[6],box[7],color);
    fillthr3D(box[6],box[7],box[5],color);
    fillthr3D(box[6],box[4],box[5],color);
    fillthr3D(box[0],box[2],box[6],color);
    fillthr3D(box[0],box[4],box[6],color);
    fillthr3D(box[1],box[5],box[7],color);
    fillthr3D(box[1],box[3],box[7],color);
    fillthr3D(box[0],box[1],box[2],color);
    fillthr3D(box[1],box[2],box[3],color);
    fillthr3D(box[0],box[1],box[4],color);
    fillthr3D(box[1],box[4],box[5],color);
    fillthr3D(box[2],box[3],box[6],color);
    fillthr3D(box[3],box[6],box[7],color);
    fillthr3D(box[6],box[7],box[4],color);
    fillthr3D(box[7],box[4],box[5],color);
    fillthr3D(box[0],box[2],box[4],color);
    fillthr3D(box[2],box[4],box[6],color);
    fillthr3D(box[1],box[5],box[3],color);
    fillthr3D(box[5],box[3],box[7],color);/*
    int color2=(((color/65535)&0xff)*2/3)<<4|(((color=256)&0xff)*2/3)<<2|((color&0xff)*2/3);
    drawline3D(box[4], box[6], color2);
    drawline3D(box[4], box[5], color2);
    drawline3D(box[5], box[7], color2);
    drawline3D(box[6], box[7], color2);
    drawline3D(box[0], box[1], color2);
    drawline3D(box[0], box[2], color2);
    drawline3D(box[0], box[4], color2);
    drawline3D(box[1], box[5], color2);
    drawline3D(box[1], box[3], color2);
    drawline3D(box[2], box[3], color2);
    drawline3D(box[2], box[6], color2);
    drawline3D(box[3], box[7], color2);*/

#else
    drawline3D(box[4], box[6], color);
    drawline3D(box[4], box[5], color);
    drawline3D(box[5], box[7], color);
    drawline3D(box[6], box[7], color);
    drawline3D(box[0], box[1], color);
    drawline3D(box[0], box[2], color);
    drawline3D(box[0], box[4], color);
    drawline3D(box[1], box[5], color);
    drawline3D(box[1], box[3], color);
    drawline3D(box[2], box[3], color);
    drawline3D(box[2], box[6], color);
    drawline3D(box[3], box[7], color);
#endif // TEST_

}


void Printc3D(int x, int y, int z, unsigned char c_dat, int color, int is)
{
    unsigned char i, j;
    c_dat -= 32;
    for (i = 0; i < 6; i++)
        for (j = 0; j < 8; j++)
            if ((Font_code[c_dat][i] << j) & 0x80)
                DrawboxM(x + i * BIG, (y + j * BIG), z, BIG/2 -0.5, rd0, color,
                         is--);

}

void Prints3D(int x, int y, int z, unsigned char *s_dat, int color)
{
    int i = 1;
    x -= strlen((const char *)s_dat) / 2 * 6 * BIG + 6 * BIG;
    x -= (strlen((const char *)s_dat) % 2) ? 3 * BIG : 0;
    while (*s_dat)
    {
        Printc3D(x += 6 * BIG, y, z - BIG / 2, *s_dat, color, i);
        s_dat++;
        i--;
    }
}
void Printc3D2(int x, int y, int z, unsigned char c_dat, int color, int is)
{
    unsigned char i, j;
    c_dat -= 32;
    for (i = 0; i < 6; i++)
        for (j = 0; j < 8; j++)
            if ((Font_code[c_dat][i] << j) & 0x80)
                DrawboxM2(x + i * BIG, (y + j * BIG), z, BIG/2 -0.5, rd1,rd0, color,
                          is--);

}

void Prints3D2(int x, int y, int z, unsigned char *s_dat, int color)
{
    int i = 1;
    x -= strlen((const char *)s_dat) / 2 * 6 * BIG + 6 * BIG;
    x -= (strlen((const char *)s_dat) % 2) ? 3 * BIG : 0;
    while (*s_dat)
    {
        Printc3D2(x += 6 * BIG, y, z - BIG / 2, *s_dat, color, i);
        s_dat++;
        i--;
    }
}

void setcam(point3D p3d, pointR3D r3d)
{
    cam0.z = p3d.z;
    camR0.rz = r3d.rz;
    cam0.y = p3d.y;
    camR0.ry = -r3d.ry;
    cam0.x = p3d.x;
    camR0.rx = -r3d.rx;
}

void seteye(point3D p3d)
{
    eye0.x = p3d.x;
    eye0.y = p3d.y;
    eye0.z = p3d.z;
}



point3D setpoint(float x, float y, float z, int c)
{
    point3D p1;
    p1.x = x;
    p1.y = y;
    p1.z = z;
    p1.c = c;
    return p1;
}

pointR3D setpointR(float x, float y, float z)
{
    pointR3D p1;
    p1.rx = x;
    p1.ry = y;
    p1.rz = z;
    return p1;
}


void setrd0(pointR3D rd)
{
    rd0=rd;
}
void setrd1(pointR3D rd)
{
    rd1=rd;
}



void drawchar(int x,int y,unsigned char ch,int c)
{
    y+=8;
    unsigned char i, j;
    ch -= 32;
    for (i = 0; i < 6; i++)
        for (j = 0; j < 8; j++)
            if ((Font_code[ch][i] << j) & 0x80)
                DrawPixel(x + i, y - j,c);
}
void drawstring(int x,int y,unsigned char* ch,int c)
{
    x-=6;
    while (*ch)
    {
        drawchar(x += 6, y, *ch, c);
        ch++;
    }

}

float D2R(float D)
{
    return 3.14159265358979323846/180*D;

}
float R2D(float D)
{
    return D/3.14159265358979323846*180;

}

int movecamR(float x,float y,float z)
{
    float iftoobig=R2D(camR0.rx);
    if(iftoobig>90&&x>0)return 0;
    if(iftoobig<-90&&x<0)return 0;
    camR0.rx+=D2R(x);
    camR0.ry+=D2R(y);
    camR0.rz+=D2R(z);
    return 1;
}
int movecam(float x,float y,float z)
{
    cam0.x+=z*sinf(camR0.ry);
    cam0.z+=z*cosf(camR0.ry);
    cam0.z-=x*sinf(camR0.ry);
    cam0.z-=y*sinf(camR0.rx);
    cam0.x+=x*cosf(camR0.ry);
    cam0.y+=y*cosf(camR0.rx);
    cam0.y+=z*sinf(camR0.rx);
    return 0;
}
