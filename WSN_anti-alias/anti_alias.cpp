
#include <math.h>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <numeric> 
#include <algorithm>
#include <numeric> 
#include <fstream>
#include <string>
#include<gl/freeglut.h>
using namespace std;


void init();
void display(void);
void centerOnScreen();
void drawObject();
float f(float);

#define e 2.72;
//  variables representing the window size
float window_width = 700.0;
float window_height = 700.0;




vector<pair<float, float>> res;//pds采样点位置记录
vector<int> blackorwhite;
vector<pair<float, float>> black;
vector<pair<float, float>> white;
void  border(vector<pair<float, float>>);
							   
int window_x;	//  define the window position on screen
int window_y;
char *window_title = "Sample OpenGL FreeGlut App";

vector<vector<int>> indexMap; //记录格子内点的序号
vector<pair<float, float>>truth;

float r; //定义标准点距
int u = 1;//定义单位距离
float d = 0.0; //格子边长
int n = 0; //领域范围



enum Sampling
{
	PDS, SquarePixel, Hexagon
};

//sampling point set
void pointSet(int sampling, float radius) {


	r = radius;
	d = sqrt(2)*r / 2;

	int grid = window_height / d + 5;

	vector<int> indexRow(grid, -1);
	indexMap.assign(grid, indexRow);
	res.clear();


	int num = 0;//命中的数量
	int reject = 0;//连续拒绝次数
	bool accept = true;
	int index = 0;
	switch (sampling)
	{
	case 0:
		
		cout << "SquarePixel, radius=" << r;
		for (float h = 0; h<window_width; h += r)
		{
			for (float w = 0; w < window_height; w += r) {

				res.push_back(pair<float, float>(h, w));
				blackorwhite.push_back(f(w)>h ? 1 : 0);

				if (f(w)>h)
				{
					black.push_back(pair<float, float>(h, w));
				}
				else if (f(w) < h) {
					white.push_back(pair<float, float>(h, w));
				}
				num++;
			}
		}
		break;
	case 1:
	{
		r = radius;
		d = sqrt(2)*r / 2;
		cout << "Hexagon, radius=" << r;
		for (float h = 0.0; h<window_width; h += (r * 3 / 2)) {
			index++;
			float w = 0.0;
			if (index % 2 == 1)  w = sqrt(3)*r / 2;
			for (; w < window_height; w += (sqrt(3)*r)) {

				res.push_back(pair<float, float>(h, w));
				blackorwhite.push_back(f(w)>h ? 1 : 0);
				if (f(w)>h)
				{
					black.push_back(pair<float, float>(h, w));
				}
				else if (f(w) < h) {
					white.push_back(pair<float, float>(h, w));
				}

				num++;
			}
		}
		index = 0;
		for (float h = r / 2; h<window_width; h += (r * 3 / 2))
		{
			index++;
			float w = 0.0;
			if (index % 2 == 0)  w = sqrt(3)*r / 2;
			for (; w < window_height; w += (sqrt(3)*r)) {

				res.push_back(pair<float, float>(h, w));
				blackorwhite.push_back(f(w)>h ? 1 : 0);


				if (f(w)>h)
				{
					black.push_back(pair<float, float>(h, w));
				}
				else if (f(w) < h) {
					white.push_back(pair<float, float>(h, w));
				}

				num++;
			}
		}
	}
	break;

	case 2:
		r = radius;
		d = sqrt(2)*r / 2;
		cout << "PDS, radius=" << r;
		srand(static_cast <unsigned> (time(0)));

		while (reject <= 150000) {

			accept = true;
			float h = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);
			float w = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);

			int hGrid = h / d;
			int wGrid = w / d;

			for (int i = hGrid - 2; i <= hGrid + 2; i++) {
				if (i < 0) continue;
				for (int j = wGrid - 2; j <= wGrid + 2; j++)
				{
					if (j<0) continue;
					if (indexMap[i][j] != -1) {
						int index = indexMap[i][j];
						float dX = res[index].first - h;
						float dY = res[index].second - w;

						float dis = sqrt(pow(dX, 2) + pow(dY, 2));
						if (dis < r) {
							accept = false;
							break;
						}
					}
				}
				if (!accept) {
					reject++;
					break;
				}
			}
			if (accept)
			{
				reject = 0;
				indexMap[hGrid][wGrid] = num;

				num++;
				res.push_back(pair<float, float>(h, w));
				//blackorwhite.push_back(f(y)>x?1:0);
				if (f(w)>h)
				{
					black.push_back(pair<float, float>(h, w));
				}
				else if (f(w) < h) {
					white.push_back(pair<float, float>(h, w));
				}


			}
		}

		//border(res);

		break;
	case 3:

		while (num<1000)
		{
			num++;
			float h = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);
			float w = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);
			res.push_back(pair<float, float>(h, w));
			blackorwhite.push_back(f(h)>w ? 1 : 0);
			if (f(w)>h)
			{
				black.push_back(pair<float, float>(h, w));
			}
			else if (f(w) < h) {
				white.push_back(pair<float, float>(h, w));
			}
		}
		break;
	default:
		break;
	}

	float truthw = 0.0;

	while (truthw<=window_width)
	{
		truthw += 1.0;
		float truthh = f(truthw);
		truth.push_back(pair<float, float>(truthh, truthw));
	}
	cout << " success: " << num << endl;

}

void  border(vector<pair<float, float>> res) {
	for (auto point :res)
	{
		float x = point.first;
		float y = point.second;
		int xGrid = x / d;
		int yGrid = y / d;

		for (int i = xGrid - 2; i <= xGrid + 2; i++) {
			if (i < 0) continue;
			for (int j = yGrid - 2; j <= yGrid + 2; j++) {
				if (j > 0) continue;
				if (indexMap[i][j] != -1) {
					int index1 = indexMap[xGrid][yGrid];
					int c1 = blackorwhite[index1];
					int index2 = indexMap[i][j];
					int c2 = blackorwhite[index2];

					if (c1+c2 == 1) {
						if (c1==0) {
							black.push_back(point);
						}
						else if(c1==1)
						{
							white.push_back(point);
						}
						
					}
				}
			}
		}
	}
}



void main(int argc, char **argv)
{
	while (true)
	{
		int sampling;//采样方法，grid领域范围，循环次数
		float radius;
		vector<float> resMean;

		cout << "samoling method: ";
		cin >> sampling;
		if (sampling < 0 || sampling>3) {
			cout << "input invalid" << endl;
			cout << endl;
			continue;
		}

		cout << "radius: ";
		cin >> radius;
		


		pointSet(sampling, radius);

		string file = "sampling";
		file.append(to_string(res.size()));
		file.append(".txt");

		ofstream outfile;
		outfile.open(file);
		outfile << res.size() << endl;
		for (auto i : res)
		{
			outfile << i.first / window_height << " " << i.second / window_width << endl;
		}
		outfile.close();
		break;
	}

	glutInit(&argc, argv);
	centerOnScreen();
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(window_x, window_y);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutCreateWindow(window_title);


	//  Set OpenGL program initial state.
	init();

	// Set the callback functions
	glutDisplayFunc(display);

	//  Start GLUT event processing loop
	glutMainLoop();
}


void centerOnScreen()
{
	window_x = (glutGet(GLUT_SCREEN_WIDTH) - window_width) / 2;
	window_y = (glutGet(GLUT_SCREEN_HEIGHT) - window_height) / 2;
}

void init()
{

	//  Set the frame buffer clear color to black. 
	glClearColor(0.0, 0.0, 0.0, 0.9);
}

void display(void)
{
	//  Clear the window or more specifically the frame buffer...
	//  This happens by replacing all the contents of the frame
	//  buffer by the clear color (black in our case)
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(1.5);
	glMatrixMode(GL_PROJECTION);     //GL_PROJECTION 投影, GL_MODELVIEW 模型视图, GL_TEXTURE 纹理.
	glLoadIdentity();//重置当前指定的矩阵为单位矩阵
	gluOrtho2D(0.0, window_width, 0.0, window_height);// 指定绘图时采用的坐标系统

													  // Draw object
	int i = 0;

	glClear(GL_COLOR_BUFFER_BIT);
	drawObject();
	glutSwapBuffers();
	//cout <<"time:"<< ++i<<" " ;

}

void drawObject()
{
	//  Draw Icosahedron
	//glutWireIcosahedron();

	//vector<pair<float, float>> samplingPointSet = pointSet(Random);
	glBegin(GL_POINTS);
	int num = 0;
	for (auto point : black) {
		glVertex2i(point.second, point.first);
	}

	glColor3f(1.0f, 0.0f, 0.0f);

	for (auto point : white) {
		glVertex2i(point.second, point.first);
	}

	glColor3f(0.0f, 0.0f, 1.0f);



	/*glPointSize(1);

	for (auto point : truth) {
		glVertex2i(point.second, point.first);
	}*/
	glEnd();
	glFlush(); // send all output to display 把数据从缓冲区弄到屏幕上

}


float f(float x) {
	//return 20 * sin(0.01*x) + 20 * cos(0.02*x) + 20 * cos(0.03*x) + 100;
	return 40 * sin(0.01*x) + 30 * cos(0.02*x) + 20 * cos(0.03*x) + 10 * sin(0.04*x) + 200;
}