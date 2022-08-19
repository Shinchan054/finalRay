//ray tracing
#include<bits/stdc++.h>
#include <windows.h>
#include <glut.h>
#include "bitmap_image.hpp"
#define pi (2*acos(0.0))
#define infinity numeric_limits<double>::infinity()
using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double length;
double radius;
double fovY;
int pixel;
int level_of_recursion;
int object_no;

ifstream myfile;
struct point
{
	double x,y,z;
};

struct ray
{
    struct point pos;
    struct point dir;
};

struct point_light
{
   point pos;
   double color[3];
};

struct spot_light
{
   point pos;
   double color[3];
   point dir;
   double cut_off;
};

vector<point_light> lights;
vector<spot_light> spotlights;

struct point pos, u ,r ,l;
struct point temp;
struct point cross_product(struct point a, struct point b)
{
    struct point res;
    res.x = a.y*b.z - a.z*b.y;
    res.y = a.z*b.x - a.x*b.z;
    res.z = a.x*b.y - a.y*b.x;
    return res;
}

struct point const_product(struct point x, double c)
{
    struct point res;
    res.x = x.x*c;
    res.y = x.y*c;
    res.z = x.z*c;
    return res;
}

struct point add(struct point a, struct point b)
{
    struct point res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}

struct point sub(struct point a, struct point b)
{
    struct point res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return res;
}

double dot_product(point a , point b)
{
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}


double radian(double deg)
{
    return (deg*pi)/180;
}


double degToAngle(double angle)
{
    return (angle*180)/pi;
}


struct point rotate(struct point a, struct point b, int dir)
{
    struct point c = cross_product(a,b);
    struct point res;
    res.x = a.x*cos(radian(3)*dir) + c.x*sin(radian(3)*dir);
    res.y = a.y*cos(radian(3)*dir) + c.y*sin(radian(3)*dir);
    res.z = a.z*cos(radian(3)*dir) + c.z*sin(radian(3)*dir);
    return res;
}


struct point normalize(struct point a)
{
    double len = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    struct point res;
    res.x = a.x/len;
    res.y = a.y/len;
    res.z = a.z/len;
    return res;
}

double vect_dist(point a,point b)
{
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

class shape
{
    public:
    struct point ref_point;
    point x1,x2,x3;
    double height, width, length,radius;
    double A,B,C,D,E,F,G,H,I,J;
    int Floorwidth,tilewidth;
    int shine;
    double color[3];
    double co_efficients[4];

    shape(){}
    virtual void draw(){}
    void setcolor(double r, double g, double b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void addAmbient(double *color, double pointColor[])
    {
        for(int i=0;i<3;i++)
        {
            color[i] += pointColor[i]*co_efficients[0];
        }
    }
    
    void addReflection(ray r,double *color, point intersect, double pointcolor[],point normal, point_light light, ray incident)
    {
       double lambert = dot_product(const_product(incident.dir, -1), normal);
       if(lambert<0)
         lambert = 0; 
       ray reflected ={intersect, sub(incident.dir, const_product(normal, 2*dot_product(incident.dir, normal)))};
       double phong = dot_product(const_product(r.dir,-1.0),reflected.dir);
       if(phong<0)
         phong = 0;
       for(int i=0;i<3;i++)
       {
         color[i]+= light.color[i]*pointcolor[i]*lambert*co_efficients[1];
         color[i]+= light.color[i]*pointcolor[i]*pow(phong,shine)*co_efficients[2]; 
       }      
    }

    void addspotReflection(ray r,double *color, point intersect, double pointcolor[],point normal, spot_light light, ray incident)
    {
        double lambert = dot_product(const_product(incident.dir, -1), normal);
       if(lambert<0)
         lambert = 0; 
       ray reflected ={intersect, sub(incident.dir, const_product(normal, 2*dot_product(incident.dir, normal)))};
       double phong = dot_product(const_product(r.dir,-1.0),reflected.dir);
       if(phong<0)
         phong = 0;
       for(int i=0;i<3;i++)
       {
         color[i]+= light.color[i]*pointcolor[i]*lambert*co_efficients[1];
         color[i]+= light.color[i]*pointcolor[i]*pow(phong,shine)*co_efficients[2]; 
       }      
    }

    void addRecursiveReflection(double *color, double pointcolor[])
    {
        for(int i=0;i<3;i++)
        {
            color[i] += pointcolor[i]*co_efficients[3];
        }
    }

    void setshine(int s){
        shine = s;
    }
    void setcoefficients(double a, double b, double c, double d){
        co_efficients[0] = a;
        co_efficients[1] = b;
        co_efficients[2] = c;
        co_efficients[3] = d;
    }
    virtual double get_intersecting_point(ray r)
    {
        return -1;
    }
    virtual double intersect(ray r, double *cur_color, int level)
    {
        return -1;
    }
    

    ~shape(){}

};

vector<shape*> shapes;


class sphere : public shape
{
    public:
    sphere(point center, double radius)
    {
        ref_point = center;
        this->radius = radius;
    }

    void draw()
    {
        glColor3f(color[0],color[1],color[2]);
        glPushMatrix();
        glTranslatef(ref_point.x,ref_point.y,ref_point.z);
        glutSolidSphere(radius,50,50);
        glPopMatrix();
    }

    double intersect(ray r, double *cur_color, int level)
    {
        double tmin;
        double a = 1.0;
        double b = 2*dot_product(sub(r.pos, ref_point), r.dir);
        double c = dot_product(sub(r.pos, ref_point), sub(r.pos, ref_point)) - radius*radius;
        double d = b*b - 4*a*c;
        if(d<0.0)
        tmin = infinity;
        else
        {
            double t1 = (-b + sqrt(d))/(2*a);
            double t2 = (-b - sqrt(d))/(2*a);
            if(t2<0.0)
            tmin = t1;
            else
            tmin = t2;
        }

        if(level==0)
        return tmin;

        point intersect_point = add(r.pos, const_product(r.dir, tmin));
        double intersect_point_color[3];
        for(int i=0;i<3;i++)
        {
            intersect_point_color[i] = this->color[i];
        }
        point normal = sub(intersect_point, ref_point);
        normal = normalize(normal);
        if(vect_dist(pos, ref_point)<=radius)
        {
           normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }

        addAmbient(cur_color, intersect_point_color);
        for(int i=0;i<lights.size();i++)
        {
            ray incident;
            incident.pos = lights[i].pos;
            incident.dir = sub(intersect_point, lights[i].pos);
            incident.dir = normalize(incident.dir);
            double t,tm = infinity;
            for(int j=0;j<shapes.size();j++)
            {
                double *dummycolor;
                t = shapes[j]->intersect(incident, dummycolor, 0);
                if(t>0.0 && t<tm)
                tm = t;
            }

            point shadowIntersect = add(incident.pos, const_product(incident.dir, tm));
            if(vect_dist(intersect_point,incident.pos)-vect_dist(shadowIntersect,incident.pos)>0.0000001)
            {
                continue;
            }
            addReflection(r, cur_color, intersect_point, intersect_point_color, normal, lights[i], incident);
        }

        if(level>=level_of_recursion)
        return tmin;

        point reflection_dir = sub(r.dir, const_product(normal, 2*dot_product(r.dir, normal)));
        reflection_dir = normalize(reflection_dir);
        ray reflected = {add(intersect_point,reflection_dir), reflection_dir};
        int nearest = INT_MAX;
        double t,tm = infinity;
        for(int j=0;j<shapes.size();j++)
        {
            double *dummycolor;
            t = shapes[j]->intersect(reflected, dummycolor, 0);
            if(t>0.0 && t<tm)
            {
                tm = t;
                nearest = j;
            }
        }

        double *reflected_color = new double[3];

        if(nearest!=INT_MAX)
        {
            tm = shapes[nearest]->intersect(reflected, reflected_color, level+1);
        }

        addRecursiveReflection(cur_color, reflected_color);

        return tmin;

    }
};

class Floor : public shape
{
    public:
    Floor(int Floorwidth, int tilewidth){
        this->Floorwidth = Floorwidth;
        this->tilewidth = tilewidth;
        ref_point.x = -Floorwidth/2.0;
        ref_point.y = -Floorwidth/2.0;
        ref_point.z = 0;
    }

    void draw()
    {
        int tiles = Floorwidth/tilewidth;
        int cnt =0;
        glBegin(GL_QUADS);
        {
            for(int i=0;i<tiles;i++)
            {
                for(int j=0;j<tiles;j++)
                {
                    glColor3f((i+j)%2,(i+j)%2,(i+j)%2);
                    glVertex3f(ref_point.x+i*tilewidth,ref_point.y+j*tilewidth,ref_point.z);
                    glVertex3f(ref_point.x+i*tilewidth,ref_point.y+j*tilewidth+tilewidth,ref_point.z);
                    glVertex3f(ref_point.x+i*tilewidth+tilewidth,ref_point.y+j*tilewidth+tilewidth,ref_point.z);
                    glVertex3f(ref_point.x+i*tilewidth+tilewidth,ref_point.y+j*tilewidth,ref_point.z);
                }
            }
        }
        glEnd();
    }

    double intersect(ray r, double *cur_color, int level)
    {
        point normal={0,0,1};
        if(dot_product(pos,normal)<=0.0)
        {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        double tmin = infinity;

        if(dot_product(r.dir,normal)!=0.0)
        {
            tmin = (-1.0)*dot_product(normal,r.pos)/dot_product(normal,r.dir);
        }
        if(tmin>0.0 && tmin<infinity)
        {
            point intersect_point = add(r.pos, const_product(r.dir, tmin));
            if(intersect_point.x < -Floorwidth/2 || intersect_point.x > Floorwidth/2 || intersect_point.y < -Floorwidth/2 || intersect_point.y > Floorwidth/2)
            {
                tmin = infinity;
            }
        }


        if(level==0)
        return tmin;

        point intersect_point = add(r.pos, const_product(r.dir, tmin));
        point ref_pos = sub(intersect_point, ref_point);
        int x = ref_pos.x/tilewidth;
        int y = ref_pos.y/tilewidth;
        double intersect_point_color[3];
        for(int i=0;i<3;i++)
        {
            intersect_point_color[i] = (x+y)%2;
        }

        addAmbient(cur_color, intersect_point_color);
        for(int i=0;i<lights.size();i++)
        {
            ray incident;
            incident.pos = lights[i].pos;
            incident.dir = sub(intersect_point, lights[i].pos);
            incident.dir = normalize(incident.dir);
            double t,tm = infinity;
            for(int j=0;j<shapes.size();j++)
            {
                double *dummycolor;
                t = shapes[j]->intersect(incident, dummycolor, 0);
                if(t>0.0 && t<tm)
                tm = t;
            }

            point shadowIntersect = add(incident.pos, const_product(incident.dir, tm));
            if(vect_dist(intersect_point,incident.pos)-vect_dist(shadowIntersect,incident.pos)>0.0000001)
            {
                continue;
            }
            addReflection(r, cur_color, intersect_point, intersect_point_color, normal, lights[i], incident);
        }

        if(level>=level_of_recursion)
        return tmin;

        point reflection_dir = sub(r.dir, const_product(normal, 2*dot_product(r.dir, normal)));
        reflection_dir = normalize(reflection_dir);
        ray reflected = {add(intersect_point,reflection_dir), reflection_dir};
        int nearest = INT_MAX;
        double t,tm = infinity;
        for(int j=0;j<shapes.size();j++)
        {
            double *dummycolor;
            t = shapes[j]->intersect(reflected, dummycolor, 0);
            if(t>0.0 && t<tm)
            {
                tm = t;
                nearest = j;
            }
        }

        double *reflected_color = new double[3];

        if(nearest!=INT_MAX)
        {
            tm = shapes[nearest]->intersect(reflected, reflected_color, level+1);
        }

        addRecursiveReflection(cur_color, reflected_color);

        return tmin;
    }
};


class triangle : public shape
{
    public:
    triangle(point p, point q, point r)
    {
        x1 = p;
        x2 = q;
        x3 = r;
    }

    void draw()
    {
        glBegin(GL_TRIANGLES);
        {
            glColor3f(color[0],color[1],color[2]);
            glVertex3f(x1.x,x1.y,x1.z);
            glVertex3f(x2.x,x2.y,x2.z);
            glVertex3f(x3.x,x3.y,x3.z);
        }
        glEnd();
    }

    double get_intersecting_point(ray r)
    {
        double eps = 0.0000001;
        point edge1,edge2,h,s,q;
        double a,f,u,v;
        edge1 = sub(x2,x1);
        edge2 = sub(x3,x1);
        h = cross_product(r.dir,edge2);
        a = dot_product(edge1,h);
        if(a > -eps && a < eps)
        return infinity;
        f = 1.0/a;
        s = sub(r.pos,x1);
        u = f*dot_product(s,h);
        if(u < 0.0 || u > 1.0)
        return infinity;
        q = cross_product(s,edge1);
        v = f*dot_product(r.dir,q);
        if(v < 0.0 || u+v > 1.0)
        return infinity;
        double t = f*dot_product(edge2,q);
        if(t > eps)
        return t;
        else
        return infinity;
    }

    double intersect(ray r, double *cur_color, int level)
    { 
        double tmin = get_intersecting_point(r);
        if(level==0)
        return tmin;

        point intersect_point = add(r.pos, const_product(r.dir, tmin));
        double intersect_point_color[3];
        for(int i=0;i<3;i++)
        {
            intersect_point_color[i] = this->color[i];
        }
        point normal = cross_product(sub(x2,x1),sub(x3,x1));
        normal = normalize(normal);
        
        if(dot_product(const_product(r.dir,-1.0),normal)<=0.0)
        {
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        addAmbient(cur_color, intersect_point_color);
        for(int i=0;i<lights.size();i++)
        {
            ray incident;
            incident.pos = lights[i].pos;
            incident.dir = sub(intersect_point, lights[i].pos);
            incident.dir = normalize(incident.dir);
            double t,tm = infinity;
            for(int j=0;j<shapes.size();j++)
            {
                double *dummycolor;
                t = shapes[j]->intersect(incident, dummycolor, 0);
                if(t>0.0 && t<tm)
                tm = t;
            }

            point shadowIntersect = add(incident.pos, const_product(incident.dir, tm));
            if(vect_dist(intersect_point,incident.pos)-vect_dist(shadowIntersect,incident.pos)>0.0000001)
            {
                continue;
            }
            addReflection(r, cur_color, intersect_point, intersect_point_color, normal, lights[i], incident);
        }

        if(level>=level_of_recursion)
        return tmin;

        point reflection_dir = sub(r.dir, const_product(normal, 2*dot_product(r.dir, normal)));
        reflection_dir = normalize(reflection_dir);
        ray reflected = {add(intersect_point,reflection_dir), reflection_dir};
        int nearest = INT_MAX;
        double t,tm = infinity;
        for(int j=0;j<shapes.size();j++)
        {
            double *dummycolor;
            t = shapes[j]->intersect(reflected, dummycolor, 0);
            if(t>0.0 && t<tm)
            {
                tm = t;
                nearest = j;
            }
        }

        double *reflected_color = new double[3];

        if(nearest!=INT_MAX)
        {
            tm = shapes[nearest]->intersect(reflected, reflected_color, level+1);
        }

        addRecursiveReflection(cur_color, reflected_color);

        return tmin;
   
    }

       
};

class general : public shape
{
    public:
    general(double a,double b,double c,double d,double e,double f,double g,double h,double i,double j,point ref_point,double height,double width,double length){
        this->A = a;
        this->B = b;
        this->C = c;
        this->D = d;
        this->E = e;
        this->F = f;
        this->G = g;
        this->H = h;
        this->I = i;
        this->J = j;
        this->ref_point = ref_point;
        this->height = height;
        this->width = width;
        this->length = length;
    }

};



void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void capture()
{
    bitmap_image image(pixel, pixel);
    
    double plane_distance = (500.0/2.0)/tan(radian(fovY/2.0));
    struct point topleft;
    struct point ll = const_product(l, plane_distance);
    struct point rr = const_product(r,500/2);
    struct point uu = const_product(u,500/2);
    struct point temp = add(ll,uu);
    topleft = add(pos,temp);
    topleft = sub(topleft,rr);
    
    double du = 500.0/pixel;
    double dv = 500.0/pixel;
    
    topleft = add(topleft, sub(const_product(r,0.5*du), const_product(u,0.5*dv)));

    for(int i=0;i<pixel;i++)
    {
        for(int j=0;j<pixel;j++)
        {
            point curpixel = add(topleft, sub(const_product(r,i*du), const_product(u,j*dv)));
            point temp = sub(curpixel, pos);
            temp = normalize(temp);
            ray r={pos,temp};
            int nearest = INT_MAX;
            int k;
            double t,tmin  = infinity;
            for(k=0;k<shapes.size();k++)
            {
                double *color = new double[3];
                t = shapes[k]->intersect(r,color,0);
                if(t>0 && t<tmin)
                {
                    tmin = t;
                    nearest = k;
                }
            }
            if(nearest!= INT_MAX)
            {
                double *color = new double[3];
                t = shapes[nearest]->intersect(r,color,1);
                for(int l=0;l<3;l++)
                {
                    if(color[l]>1.0)color[l]=1;
                    if(color[l]<0.0)color[l]=0;
                }
            image.set_pixel(i,j,(int)round(color[0]*255.0),(int)round(color[1]*255),(int)round(color[2]*255));
            }

            }
        }

         cout<<"saving image"<<endl;
    image.save_image("image.bmp");
    image.clear();
    }


void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
		  l = rotate(l,u,-1);
          r = rotate(r,u,-1);	
			break;
        case '2':
         l = rotate(l,u,1);
         r =  rotate(r,u,1);
         break;
        case '3':
        l = rotate(l,r,-1);
        u = rotate(u,r,-1);
        break;
        case '4':
        l = rotate(l,r,1);
        u = rotate(u,r,1);
        break;
        case '5':
        u= rotate(u,l,-1);
        r= rotate(r,l,-1);
        break;
        case '6':
        u = rotate(u,l,1);
        r = rotate(r,l,1);
        break;
        case '0':
        capture();
        break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			temp = const_product(l,2);
            pos = sub(pos,temp);
			break;
		case GLUT_KEY_UP:		// up arrow key
			temp = const_product(l,2);
            pos = add(pos,temp);
			break;

		case GLUT_KEY_RIGHT:
			temp = const_product(r,2);
            pos = add(pos,temp);
			break;
		case GLUT_KEY_LEFT:
			temp = const_product(r,2);
            pos = sub(pos,temp);
			break;

		case GLUT_KEY_PAGE_UP:
            temp = const_product(u,2);
            pos = add(pos,temp);
			break;
		case GLUT_KEY_PAGE_DOWN:
            temp = const_product(u,2);
            pos = sub(pos,temp);
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
        if(length>0)
          length--;
        if(radius<30)
        radius++;  
			break;
		case GLUT_KEY_END:
        if(length<30)
            length++;
        if(radius>0)
            radius--;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}




void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x,pos.y,pos.z, pos.x+l.x,pos.y+l.y,pos.z+l.z, u.x,u.y,u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

     for(int i=0;i<shapes.size();i++)
     {
        shapes[i]->draw();       
     }

     for(int i=0;i<lights.size();i++)
     {
        
            glPushMatrix();
            glColor3f(lights[i].color[0],lights[i].color[1],lights[i].color[2]);
            glTranslatef(lights[i].pos.x,lights[i].pos.y,lights[i].pos.z);
            glutSolidSphere(2,100,100);
            glPopMatrix();    
     }

     for(int i=0;i<spotlights.size();i++)
     {
            glPushMatrix();
            glColor3f(spotlights[i].color[0],spotlights[i].color[1],spotlights[i].color[2]);
            glTranslatef(spotlights[i].pos.x,spotlights[i].pos.y,spotlights[i].pos.z);
            glutSolidSphere(5,100,100);
            glPopMatrix();    
     }

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void loadData(){
        myfile>>level_of_recursion;
        myfile>>pixel;
        myfile>>object_no;
        cout<<object_no<<endl;
        for(int i=0;i<object_no;i++)
        {
            string str;
            
            myfile>>str;
            if(str=="sphere")
            {
              shape *temp;
               point center;
               double radius;
               double color[3];
               myfile>>center.x>>center.y>>center.z;
               myfile>>radius;
                temp= new sphere(center, radius);
                myfile>>color[0]>>color[1]>>color[2];
                temp->setcolor(color[0],color[1],color[2]);
                double coefficients[4];
                myfile>>coefficients[0]>>coefficients[1]>>coefficients[2]>>coefficients[3];
                temp->setcoefficients(coefficients[0],coefficients[1],coefficients[2],coefficients[3]);
                double shine;
                myfile>>shine;
                temp->setshine(shine);
                shapes.push_back(temp);
            }
            else if(str=="triangle")
            {
                shape *temp;
                point p1,p2,p3;
                double color[3];
                myfile>>p1.x>>p1.y>>p1.z;
                myfile>>p2.x>>p2.y>>p2.z;
                myfile>>p3.x>>p3.y>>p3.z;
                temp= new triangle(p1,p2,p3);
                myfile>>color[0]>>color[1]>>color[2];
                temp->setcolor(color[0],color[1],color[2]);
                double coefficients[4];
                myfile>>coefficients[0]>>coefficients[1]>>coefficients[2]>>coefficients[3];
                temp->setcoefficients(coefficients[0],coefficients[1],coefficients[2],coefficients[3]);
                double shine;
                myfile>>shine;
                temp->setshine(shine);
                shapes.push_back(temp);

            }
            else if(str=="general")
            {
                shape *temp;
                double a,b,c,d,e,f,g,h,i,j;
                point ref;
                double height,width,length;
                double color[3];
                double coefficients[4];
                myfile>>a>>b>>c>>d>>e>>f>>g>>h>>i>>j;
                myfile>>ref.x>>ref.y>>ref.z>>length>>width>>height;
                temp= new general(a,b,c,d,e,f,g,h,i,j,ref,height,width,length);
                myfile>>color[0]>>color[1]>>color[2];
                temp->setcolor(color[0],color[1],color[2]);
                myfile>>coefficients[0]>>coefficients[1]>>coefficients[2]>>coefficients[3];
                temp->setcoefficients(coefficients[0],coefficients[1],coefficients[2],coefficients[3]);
                double shine;
                myfile>>shine;
                temp->setshine(shine);
                //shapes.push_back(temp);
                
            }
        }

        shape *temp;
        temp = new Floor(1000,20);
        temp->setcoefficients(0.4,0.2,0.2,0.2);
        temp->setshine(10);
        shapes.push_back(temp);
        
        cout<<shapes.size()<<endl;
        int point;
        myfile>>point;
        for(int i=0;i<point;i++)
        {
            point_light p;
            myfile>>p.pos.x>>p.pos.y>>p.pos.z;
            myfile>>p.color[0]>>p.color[1]>>p.color[2];
            lights.push_back(p);
        }
      
        cout<<lights.size()<<endl;
        int spot;
        myfile>>spot;
        for(int i=0;i<spot;i++)
        {
            spot_light s;
            myfile>>s.pos.x>>s.pos.y>>s.pos.z;
            myfile>>s.color[0]>>s.color[1]>>s.color[2];
            myfile>>s.dir.x>>s.dir.y>>s.dir.z;
            myfile>>s.cut_off;
            s.dir = normalize(s.dir);
            spotlights.push_back(s);
        }
    }



void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
    fovY = 80;
    length = 20.0;
    radius = 10.0;
    pos.x=100;
    pos.y=100;
    pos.z=0;
    u.x = 0;
    u.y=0;
    u.z=1;
    r.x = -1/sqrt(2);
    r.y = 1/sqrt(2);
    r.z = 0;
    l.x = -1/sqrt(2);
    l.y = -1/sqrt(2);
    l.z = 0;
    myfile.open("scene.txt",ios::in);
	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();
   
	//give PERSPECTIVE parameters
	gluPerspective(fovY,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();
    loadData();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
