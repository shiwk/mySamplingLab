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
#include <fstream>  
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
double angle(pair<float, float>);



//  define the window position on screen
int window_x;
int window_y;

//  variables representing the window size
float window_width =350.0;
float window_height = 350.0;
int N = 0;
int numRandom = 0;
int rejectTime = 0;
// 剩余点数，以及连续没有消点成功次数
int reject = 0;


//  variable representing the window title
char *window_title = "Sample OpenGL FreeGlut App";


float ratio = 0.75;
float R = 0.0;//定义标准点距
																		 //float u = 1.0;//定义单位距离
float d =0.0; //格子边长

//grid num
int w = 0;
int h =0;

//Random中用来记录每个点的grid位置
vector < pair<float, float>> points ;
vector<vector<pair<float, float>>> rows;
vector<vector<vector<pair<float, float>>>> samplingRandomPointSet;
vector<vector<vector<pair<float, float>>>> initialPointSet;

int deletePoints = 0;


//sampling point set
void pointSet() {

	srand(static_cast <unsigned> (time(0)));

	for (int i = 0; i < numRandom; i++) {
		float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);
		float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);

		int xGrid = x / d;
		int yGrid = y / d;

		initialPointSet[xGrid][yGrid].push_back(pair<float, float>(x, y));
	}
	


}

//-------------------------------------------------------------------------
//  Program Main method.
//-------------------------------------------------------------------------
void main(int argc, char **argv)
{
	cout << "numRandom:";
	cin >> numRandom;
	cout << "N:";
	cin >> N;

	/*cout << "rejectTime:";
	cin >> rejectTime;*/

	float D = 2 * sqrt(window_width*window_height / (2 * N*sqrt(3)));
	R = ratio * D;
	d = sqrt(2)*R / 2;

	float min = D;

	w = static_cast <int> (window_width / d) + 5;
	h = static_cast <int> (window_height / d) + 5;

	points = {};
	rows = vector<vector<pair<float, float>>>(w, points);
	initialPointSet = vector<vector<vector<pair<float, float>>>>(h, rows);


	pointSet();
	samplingRandomPointSet = initialPointSet;

	cout << "R:" << R << "  time:" << currentDateTime() << endl;

	int numPoints = 0;

	int i = 0;//投掷点数



	numPoints = 0;
	deletePoints = 0;
	samplingRandomPointSet = initialPointSet;
	reject = 0;
	R = D*ratio;
	cout << "ratio:" << ratio << "  R:" << R << endl;


	int tmp = deletePoints;

	//samplingRandomPointSet = initialPointSet;
	while (true) {
		i++;
		deletePoint();

		//if (deletePoints == tmp) {
		//	reject++;
		//	if (reject == rejectTime) {
		//		cout << currentDateTime() << endl;
		//		//system("pause");
		//		cout << "deletePoints: " << deletePoints << endl;
		//		break;
		//	}
		//	continue;
		//}
		//else {
		//	reject = 0;
		//	//cout << i << " " << deletePoints << endl;
		//	tmp = deletePoints;
		//}

		if (deletePoints >= numRandom - N-100) {
			cout << currentDateTime() << endl;
			//system("pause");
			break;
		}
		//else cout << deletePoints << endl;

	}
	numPoints = numRandom - deletePoints;
	string file = "sampling";
	file.append(to_string(numPoints));
	file.append(".txt");

	ofstream outfile;
	outfile.open(file);
	outfile << numPoints << endl;
	for (auto i : samplingRandomPointSet) {
		for (auto j : i)
		{
			for (auto point : j) {
				outfile << point.first / window_height << " " << point.second / window_width << endl;


				float x = point.first;
				float y = point.second;
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
								if (dX == 0.0 && dY == 0.0) continue;
								float dis = sqrt(pow(dX, 2) + pow(dY, 2));

								min = min < dis ? min : dis;
							}
						}
					}
				}

			}

		}
	}
	outfile.close();

	cout << "min:" << min << endl;

	float realR = 2 * sqrt(window_width*window_height / (2 * (numPoints)*sqrt(3)));
	cout << "realR:" << realR << endl;
	cout << "relative radius:" << min / realR << endl;


	cout << "point num:" << numPoints << endl << endl;
	ratio = (static_cast <float>(D - realR) / D + 1)*ratio;
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


void init()
{

	//  Set the frame buffer clear color to black. 
	glClearColor(1.0, 1.0, 1.0, 0.9);
}




void display(void)
{
	//  Clear the window or more specifically the frame buffer...
	//  This happens by replacing all the contents of the frame
	//  buffer by the clear color (black in our case)
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.0f, 0.0f);
	glPointSize(2.5);
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
	for (auto i : samplingRandomPointSet) {
		for (auto j : i)
		{
			for (auto point : j) {
				glVertex2i(point.second, point.first);
				num++;
			}

		}
	}

	cout << "point num:" << num;

	
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
						int s = floor(a / (M_PI /3 ));

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
		}
	}
	//cout << " Neighbor:" << numNei << endl;
	if(numNei==1) samplingRandomPointSet[xGrid][yGrid].push_back(pair<float, float>(x, y));
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

				deletePoints += (i.size() - 1);
			}
		}
	}
	
}



//计算角度，即分区
double angle(pair<float, float> p) {
	double dX = p.first;
	double dY = p.second;

	float vectorA[] = { -1,0 };
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