#include <gl/freeglut.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <cstdlib>
#include <Windows.h>
#include <algorithm>
#include <numeric> 
#include <string>
using namespace std;




// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}


#define e 2.72;

void init();
void display(void);
void centerOnScreen();
void drawObject();
void deletePoint();
double angle(pair<float,float>);
vector<pair<float, float>> nei(float, float);
pair<float, float> xAndy(float a1, float b1, float c1, float a2, float b2, float c2);
void maximum(float,float);
bool ban(float,float);


//  define the window position on screen
int window_x;
int window_y;

//  variables representing the window size
float window_width = 350.0;
float window_height = 350.0;
int randomNum = 300000;
int N = 1000;
int maxN = 100000;
int maxTime = 0;
// 剩余点数，以及连续没有消点成功次数
int remainingNum = 0;
int reject = 0;
int  deletePoints = 0;
int numPoint = 0;

//  variable representing the window title
char *window_title = "Sample OpenGL FreeGlut App";


float ratio = 0.78;
float R = ratio*2*sqrt(window_width*window_height/(2*N*sqrt(3)));//定义标准点距
//float u = 1.0;//定义单位距离
float d = sqrt(2)*R / 2; //格子边长


int w = static_cast <int> (window_width / d)+5;
int h = static_cast <int> (window_height / d) + 5;

//Random中用来记录每个点的grid位置
vector < pair<float,float>> points = {};
vector<vector<pair<float, float>>> rows (w,points);
vector<vector<vector<pair<float, float>>>> samplingRandomPointSet(h, rows);

vector <bool> flags (w,true);
vector<vector<bool>> gridFlag(h,flags);
vector<vector<bool>> gridAdded(h, flags);
//sampling point set
void pointSet() {

	srand(static_cast <unsigned> (time(0)));
			
	for (int i = 0; i < randomNum; i++) {
		float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);
		float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);

		int xGrid = x / d;
		int yGrid = y / d;

		samplingRandomPointSet[xGrid][yGrid].push_back(pair<float, float>(x, y));
		gridFlag[xGrid][yGrid] = false;
	}
	
}

//-------------------------------------------------------------------------
//  Program Main method.
//-------------------------------------------------------------------------
void main(int argc, char **argv)
{
	cout << "R:" << R<<"  time:" << currentDateTime() <<endl;
	pointSet();

	int i = 0;
	int tmp = deletePoints;
	int rejectTime = 100000;

	while (true) {
		i++;
		deletePoint();

		if (deletePoints == tmp) {
			reject++;
			if (reject == rejectTime) {
				cout << currentDateTime() << endl;
				//system("pause");
				cout << "deletePoints: " << deletePoints << endl;
				break;
			}
			continue;
		}
		else {
			reject = 0;
			//cout << i << " " << deletePoints << endl;
			tmp = deletePoints;
		}
	}
	numPoint = randomNum - deletePoints;
	cout << numPoint<<endl;
	
	
	//  Connect to the windowing system + create a window
	//  with the specified dimensions and position
	//  + set the display mode + specify the window title.
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

//-------------------------------------------------------------------------
//  Set OpenGL program initial state.
//-------------------------------------------------------------------------
void init()
{
	//  Set the frame buffer clear color to black. 
	glClearColor(1.0, 1.0, 1.0, 0.9);
}

//-------------------------------------------------------------------------
//  This function is passed to glutDisplayFunc in order to display 
//  OpenGL contents on the window.
//-------------------------------------------------------------------------
void display(void)
{
	//  Clear the window or more specifically the frame buffer...
	//  This happens by replacing all the contents of the frame
	//  buffer by the clear color (black in our case)
	glClear(GL_COLOR_BUFFER_BIT);
	//glColor3f(0.0f, 0.0f, 0.0f);
	glPointSize(1.5);
	glMatrixMode(GL_PROJECTION);     //GL_PROJECTION 投影, GL_MODELVIEW 模型视图, GL_TEXTURE 纹理.
	glLoadIdentity();//重置当前指定的矩阵为单位矩阵
	gluOrtho2D(0.0, window_width, 0.0, window_height);// 指定绘图时采用的坐标系统
							
													  // Draw object
	int i = 0;

	/*while (true)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		drawObject();
		glutSwapBuffers();
		maximum();
	}*/
	glClear(GL_COLOR_BUFFER_BIT);
	drawObject();
	glutSwapBuffers();
	//cout <<"time:"<< ++i<<" " ;

	
		
}

//-------------------------------------------------------------------------
//  Draws our object.
//-------------------------------------------------------------------------
void drawObject()
{
	//  Draw Icosahedron
	//glutWireIcosahedron();

	//vector<pair<float, float>> samplingPointSet = pointSet(Random);
	int num = 0;
	glBegin(GL_POINTS);
	/*for (auto i : samplingRandomPointSet) {
		for (auto j: i )
		{
			for (auto point : j) {
				glVertex2i(point.second, point.first);
				num++;
			}
		}
	}*/
	for (int  i = 0; i < samplingRandomPointSet.size(); i++)
	{
		for (int j = 0; j < samplingRandomPointSet[i].size(); j++) {
			for (auto point: samplingRandomPointSet[i][j])
			{
				maximum(point.second, point.first);

				if (!gridAdded[i][j]) {
					glColor3f(1.0f, 0.0f, 0.0f);
					glPointSize(3.0);
				}
				else {
					glColor3f(0.0f, 0.0f, 0.0f);
					glPointSize(1.5);

				}
				glVertex2i(point.second, point.first);
				num++;
			}
		}
	}
	cout <<numPoint<<"before,"<< num << " after" << endl;
	glEnd();
	glFlush(); // send all output to display 把数据从缓冲区弄到屏幕上

}

//-------------------------------------------------------------------------
//  This function sets the window x and y coordinates
//  such that the window becomes centered
//-------------------------------------------------------------------------
void centerOnScreen()
{
	window_x = (glutGet(GLUT_SCREEN_WIDTH) - window_width) / 2;
	window_y = (glutGet(GLUT_SCREEN_HEIGHT) - window_height) / 2;
}



void deletePoint() {
	float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);
	float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);

	int xGrid = x / d;
	int yGrid = y / d;

	vector<pair<float, float>> part = {};
	vector<vector<pair<float, float>>> parts(6, part);

	int numNei = 0;
	for (int i = xGrid - 2; i <= xGrid + 2; i++) {
		if (i < 0 || i>h) continue;
		for (int j = yGrid - 2; j <= yGrid + 2; j++)
		{
			if (j < 0 || j>w) continue;
			if (samplingRandomPointSet[i][j].empty()) continue;
			else {
				for (auto it = samplingRandomPointSet[i][j].begin(); it != samplingRandomPointSet[i][j].end(); it++)
				{
					float dX = (*it).first - x;
					float dY = (*it).second - y;

					float dis = sqrt(pow(dX, 2) + pow(dY, 2));

					if (dis < R) {
						numNei++;
						pair<float, float> p(dX, dY);
						double a = angle(p);
						int s = floor(a / (M_PI / 3));

						if (s < 0 || s>5) {
							cout << "error:p=" << s << "  angle:" << a << endl;
							system("pause");
						}

						parts[s].push_back(*it);

						it = samplingRandomPointSet[i][j].erase(it);


						if (it == samplingRandomPointSet[i][j].end()) {
							break;
						}

					}
				}
			}
			if(samplingRandomPointSet[i][j].empty())
				gridFlag[i][j] = true;
		}

	}


	//cout << " Neighbor:"<<numNei << endl;
	if (numNei == 1) {
		samplingRandomPointSet[xGrid][yGrid].push_back(pair<float, float>(x, y));
		gridFlag[xGrid][yGrid] = false;
	}
	else
	{
		for (auto i : parts)
		{
			float firstSum = 0.0, secondSum = 0.0;
			for (auto j : i) {
				firstSum += j.first;
				secondSum += j.second;
			}
			if (firstSum != 0.0 || secondSum != 0.0) {
				float hPos = firstSum / i.size();
				float wPos = secondSum / i.size();

				int hGrid = hPos / d;
				int wGrid = wPos / d;

				samplingRandomPointSet[hGrid][wGrid].push_back(pair<float, float>(hPos, wPos));
				gridFlag[hGrid][wGrid] = false;
				deletePoints += (i.size() - 1);

			}
		}
	}
}



//计算角度，即分区
double angle(pair<float, float> p) {
	double dX = p.first;
	double dY = p.second;

	float vectorA[] = {-1,0 };
	float vectorB[] = { dX, dY };
	double dotProduct = inner_product(vectorA, vectorA + 2, vectorB, 0.0);

	double cos = dotProduct / (1 * sqrt(pow(dX, 2) + pow(dY, 2)));

	double res = acos(cos);
	if (dY >= 0) {
		//cout << "angle:" << res << endl;
		return res;
	}
	else {
		//cout << "angle:" << res << endl;
		return 2 * M_PI - res;
	}

}

vector<pair<float, float>> nei(float x, float y) {

	int xGrid = x / d;
	int yGrid = y / d;

	int neis = 0;
	vector<double> angles = {};
	vector<pair<float, float>> res = {};

	for (int i = xGrid - 3; i <= xGrid + 3; i++) {
		if (i < 0 || i>h) continue;
		for (int j = yGrid - 3; j <= yGrid + 3; j++)
		{
			if (j < 0 || j>w) continue;
			//if (samplingRandomPointSet[i][j].empty()) continue;

			else {
				for (auto it = samplingRandomPointSet[i][j].begin(); it != samplingRandomPointSet[i][j].end(); it++)
				{
					float dX = (*it).first - x;
					float dY = (*it).second - y;

					float dis = sqrt(pow(dX, 2) + pow(dY, 2));

					if (dis > R && dis < 2 * R) {

						neis++;
						pair<float, float> p((*it).first, (*it).second);
						pair<float, float> d(dX, dY);
						double a = angle(d);
						vector<pair<float, float>>::iterator itPoint = res.begin();
						vector<double>::iterator itAngle = angles.begin();

						while (true)
						{
							if (itAngle == angles.end()) {
								angles.push_back(a);
								res.push_back(p);
								break;
							}
							else
							{
								if (*itAngle < a) {
									itAngle++;
									itPoint++;
									continue;
								}
								else
								{
									res.insert(itPoint, p);
									angles.insert(itAngle, a);
									break;
								}

							}

						}

					}
				}
			}
		}
	}
	return res;
}


void maximum(float x,float y) {
	/*float x, y;
	while (true)
	{
		srand(static_cast <unsigned> (time(0)));

		float xGrid = static_cast <int> (rand()) / static_cast <int> (RAND_MAX / h);
		float yGrid = static_cast <int> (rand()) / static_cast <int> (RAND_MAX / w);
		if (xGrid < 0 || xGrid >= h || yGrid < 0 || yGrid >= w) continue;

		if (!samplingRandomPointSet[xGrid][yGrid].empty()) {
			x = samplingRandomPointSet[xGrid][yGrid][0].first;
			y = samplingRandomPointSet[xGrid][yGrid][0].second;
			break;
		}
	}*/
	
	vector<pair<float, float>> neighbors=nei(x,y);
	//cout << " neighbors:"<<neighbors.size() << endl;


	vector<pair<float, float>> ::const_iterator it1=neighbors.begin();
	vector<pair<float, float>> ::const_iterator it2 = neighbors.end()-1;
	float x1 = x, y1 = y;
	

	if (neighbors.size() == 0) { 
		cout << "no neighbors";
		return;
	}
	if (neighbors.size() == 1) {
		cout << "1 neighbors";
		return;
	}

	while (it1 != neighbors.end())
	{
		if (it2 == neighbors.end()) it2 = neighbors.begin();
		
		float x2 = (*it1).first;
		float y2 = (*it1).second;
		float x3 = (*it2).first;
		float y3 = (*it2).second;
		float a1 = x2 - x1,b1 = y2 - y1, c1 = -((x1+x2)*(x2-x1)/2+(y2-y1)*(y1+y2)/2);
		float a2 = x3 - x1, b2 = y3 - y1, c2 = -((x1 + x3)*(x3 - x1) / 2 + (y3 - y1)*(y1 + y3) / 2);

		it1++;
		it2++;

		if (a1*b2 == a2*b1) {
			cout << x << " "<<y << " "<<endl;
			for (auto n : neighbors) cout << n.first << "," << n.second << "  ";
			
			cout <<endl<< x2 << " " << y2 << " " << x3 << " " << y3;
			cout << " line ";
			system("pause");

			continue;
		}

		pair<float, float> excenter = xAndy(a1,b1,c1,a2,b2,c2);

		//cout <<"(" <<excenter.first << "," << excenter.second <<")   ";
		float dX = excenter.first - x;
		float dY = excenter.second - y;

		float dis = sqrt(pow(dX, 2) + pow(dY, 2));
		if (dis > R) {
			int xGrid = excenter.first / d;
			int yGrid = excenter.second / d;
			if (excenter.first < 0 || excenter.first> window_height || excenter.second < 0 || excenter.second > window_width) continue;
			if (gridFlag[xGrid][yGrid]) {
				
				if (ban(excenter.first, excenter.second)) break;

				samplingRandomPointSet[xGrid][yGrid].push_back(excenter);
				gridFlag[xGrid][yGrid] = false;
				gridAdded[xGrid][yGrid] = false;
				//cout << "maxTime:" << maxTime << endl;
				cout << "("<<x << "," << y << ")" ;
				//for (auto n : neighbors) cout << "(" << n.first << "," << n.second << ")" <<endl;
				cout <<endl<< "add(" << excenter.first << "," << excenter.second << ")" << endl<<endl;
				//cout << "(" << (*it1).first << "," << (*it1).second << ") " << "(" << (*it2).first << "," << (*it2).second << ")" << endl;
				//system("pause");

			}
		}
		
	}
}


bool ban(float x,float y) {
	bool ban = false;
	int xGrid = x / d;
	int yGrid = y / d;

	for (int i = xGrid - 2; i <= xGrid + 2; i++) {
		if (i < 0 || i>h) continue;
		for (int j = yGrid - 2; j <= yGrid + 2; j++)
		{
			if (j < 0 || j>w) continue;
			if (samplingRandomPointSet[i][j].empty()) continue;
			else {
				for (auto it = samplingRandomPointSet[i][j].begin(); it != samplingRandomPointSet[i][j].end(); it++)
				{
					float dX = (*it).first - x;
					float dY = (*it).second - y;

					float dis = sqrt(pow(dX, 2) + pow(dY, 2));

					if (dis < R) {
						return true;
					}
				}
			}
		}
	}

	return false;
}

pair<float, float> xAndy(float a1, float b1, float c1, float a2, float b2, float c2) {
	pair<float, float> res;
	float x = 0.0, y = 0.0;
	if (a1 == 0) {
		if (b1 != 0) {
			y = -c1 / b1;
			if (a2 != 0) {
				x = (b2*y - c2) / a2;
			}
		}
	}
	else if (a2 == 0) {
		if (b2 != 0) {
			y = -c2 / b2;
			if (a1 != 0) {
				x = (b1*y - c1) / a1;
			}
		}
	}
	else {
		y = (a1*c2 - a2*c1) / (a2*b1 - a1*b2);
		x = -(a2*b1*y + a2*c1) / (a1*a2);
	}


	res.first = x;
	res.second = y;


	if (res.first == 0.0&&res.second == 0.0) {
		cout << "error";
		system("pause");
	}
	return res;
}