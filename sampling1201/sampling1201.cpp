// sampling1201.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <ctime>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

struct Point
{
	float x;
	float y;
	float z;
	Point(float w, float h, float l) :x(w), y(h), z(l) {};
};

vector<Point> vertexUp(Point);
vector<Point> vertexDown(Point);
vector<Point> centers(vector<Point>, vector<Point>);
float distance(Point , Point );
int partition(Point, vector<Point>);
void deletePoint();


//  variables representing the window size
float window_width = 350.0;
float window_height = 350.0;
float window_length = 350.0;

int N = 1000;
int numRandom = 300000;
int rejectTime = 0;
// 剩余点数，以及连续没有消点成功次数
int reject = 0;



float ratio = 0.75;
float R = 0.0;//定义标准点距
			  //float u = 1.0;//定义单位距离
float d = 0.0; //格子边长

//grid num
int w = 0;
int h = 0;
int l = 0;

int deletePoints = 0;

//Random中用来记录每个点的grid位置
vector < Point> points;
vector<vector<Point>> rows;
vector<vector<vector<Point>>> squares;
vector<vector<vector<vector<Point>>>> samplingRandomPointSet;




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



//sampling point set
void pointSet() {

	srand(static_cast <unsigned> (time(0)));

	for (int i = 0; i < numRandom; i++) {
		float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);
		float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);
		float z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_length);


		int xGrid = x / d;
		int yGrid = y / d;
		int zGrid = z / d;

		samplingRandomPointSet[xGrid][yGrid][zGrid].push_back(Point(x,y,z));

	}
	cout << "Random Ok" << endl;
}

int main()
{
	
	cout << "rejectTime:";
	cin >> rejectTime;

	float D = 2 * pow(window_width*window_height*window_length / (4 * N*sqrt(2)),1.0/3);
	R = ratio * D;
	d = sqrt(3)*R / 3;


	float min = D;

	w = static_cast <int> (window_width / d) + 5;
	h = static_cast <int> (window_height / d) + 5;
	l = static_cast <int> (window_length / d) + 5;


	//initialize sampling points
	points = {};
	rows = vector<vector<Point>>(w, points);
	squares = vector<vector<vector<Point>>>(h, rows);
	samplingRandomPointSet = vector<vector<vector<vector<Point>>>>(l, squares);

	pointSet();
	
	cout << "R:" << R << "  time:" << currentDateTime() << endl;

	int i = 0;
	int tmp = deletePoints;
	while (true) {
		i++;

		deletePoint();

		if (deletePoints == tmp) {
			reject++;
			if (reject == rejectTime) {
				cout << "R:" << R << "  time:" << currentDateTime() << endl;
				//system("pause");
				cout << i << " " << deletePoints << endl;
				break;
			}
			continue;
		}
		else {

			reject = 0;
			//cout << i << " " << deletePoints << endl;
			tmp = deletePoints;
		}
		cout<<"delete "<<deletePoints<<endl;

	}


    return 0;
}

void deletePoint() {

	srand(static_cast <unsigned> (time(0)));

	float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_height);
	float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_width);
	float z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / window_length);

	Point p(x, y, z);

	int xGrid = x / d;
	int yGrid = y / d;
	int zGrid = z / d;


	vector<Point> part = {};
	vector<vector<Point>> parts(20, part);

	vector<Point> up = vertexUp(p);
	vector<Point> down = vertexDown(p);

	/*for (auto i : up) cout << "(" << i.x<<","<<i.y<<","<<i.z<< ")";
	cout << endl;
	for (auto i : down) cout << "(" << i.x << "," << i.y << "," << i.z << ")";
	cout << endl;*/

	vector<Point> center = centers(up, down);

	//for (auto i : center) cout << "(" << i.x << "," << i.y << "," << i.z << ")";



	int numNei = 0;
	for (int i = xGrid - 2; i <= xGrid + 2; i++) {
		if (i < 0 || i>h) continue;
		for (int j = yGrid - 2; j <= yGrid + 2; j++)
		{
			if (j < 0 || j>w) continue;
			for (int t = zGrid - 2; t < zGrid + 2; t++)
			{
				if (t < 0 || t>l) continue;

				if (samplingRandomPointSet[i][j][t].empty()) continue;
				else {
					for (auto it = samplingRandomPointSet[i][j][t].begin(); it != samplingRandomPointSet[i][j][t].end(); it++)
					{
						float dX = (*it).x - x;
						float dY = (*it).y - y;
						float dZ = (*it).z - z;

						float dis = sqrt(pow(dX, 2) + pow(dY, 2)+pow(dZ,2));

						if (dis < R) {
							numNei++;
						
							int s = partition(*it,center);

							if (s < 0 || s>20) {
								cout << "error:partition=" << s << endl;
								cout << "point:(" << (*it).x << "," << (*it).y << "," << (*it).z << ")" << endl;
								for (auto i : up) cout << "(" << i.x<<","<<i.y<<","<<i.z<< ")";
								cout << endl;
								for (auto i : down) cout << "(" << i.x << "," << i.y << "," << i.z << ")";
								cout << endl;
								for (auto i : center) cout << "(" << i.x << "," << i.y << "," << i.z << ")";
								cout << endl;

								system("pause");
							}

							parts[s].push_back(*it);

							it = samplingRandomPointSet[i][j][t].erase(it);


							if (it == samplingRandomPointSet[i][j][t].end()) {
								break;
							}

						}
					}
				}
			}
		}
	}
	//cout << " Neighbor:" << numNei << endl;
	if (numNei == 1) samplingRandomPointSet[xGrid][yGrid][zGrid].push_back(p);
	else
	{
		for (auto i : parts)
		{
			float firstSum = 0.0, secondSum = 0.0, thirdSum=0.0;
			for (auto j : i) {
				firstSum += j.x;
				secondSum += j.y;
				thirdSum += j.z;
			}
			if (firstSum != 0.0 || secondSum != 0.0) {
				float newX = firstSum / i.size();
				float newY = secondSum / i.size();
				float newZ = thirdSum / i.size();

				int xGrid = newX / d;
				int yGrid = newY / d;
				int zGrid = newZ / d;

				Point newPoint(newX,newY,newZ);
				samplingRandomPointSet[xGrid][yGrid][zGrid].push_back(newPoint);

				deletePoints += (i.size() - 1);
			}
		}
	}

}

vector<Point>  vertexUp(Point p) {

	//12 vertexs
	vector<Point> vertexUp;

	float x = p.x;
	float y = p.y;
	float z = p.z;

	float r = R*sqrt(3) / 2;
	vertexUp.push_back(Point(x,y,z+R));
	vertexUp.push_back(Point(x,y+r, z + R/2));
	vertexUp.push_back(Point(x+r*cos(M_PI/10), y + r*sin(M_PI/10), z + R / 2));
	vertexUp.push_back(Point(x + r*cos(M_PI*3 / 10), y - r*sin(M_PI * 3 / 10), z + R / 2));
	vertexUp.push_back(Point(x - r*cos(M_PI * 3 / 10), y - r*sin(M_PI * 3 / 10), z + R / 2));
	vertexUp.push_back(Point(x - r*cos(M_PI / 10), y + r*sin(M_PI * 1 / 10), z + R / 2));

	return vertexUp;

}

vector<Point>  vertexDown(Point p) {

	//12 vertexs
	vector<Point> vertexDown;

	float x = p.x;
	float y = p.y;
	float z = p.z;

	float r = R*sqrt(3) / 2;


	vertexDown.push_back(Point(x, y, z - R));
	vertexDown.push_back(Point(x + r*cos(M_PI * 3 / 10), y + r*sin(M_PI * 3 / 10), z - R / 2));
	vertexDown.push_back(Point(x + r*cos(M_PI / 10), y - r*sin(M_PI / 10), z - R / 2));
	vertexDown.push_back(Point(x, y - r, z - R / 2));
	vertexDown.push_back(Point(x - r*cos(M_PI / 10), y - r*sin(M_PI / 10), z - R / 2));
	vertexDown.push_back(Point(x - r*cos(M_PI * 3 / 10), y + r*sin(M_PI * 3 / 10), z - R / 2));

	return vertexDown;
}


int partition(Point p,vector<Point>center ) {
	float min = R;
	int res = -1;
	for (int i = 0; i < center.size(); i++)
	{
		float dis = distance(p, center[i]);
		if ( dis< min) {
			min = dis;
			res = i;
		}
	}
	return res;
}

float distance(Point p1, Point p2) {
	float dX = p1.x - p2.x;
	float dY = p1.y - p2.y;
	float dZ = p1.z - p2.z;

	float dis = sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2));
	return dis;
}


vector<Point> centers(vector<Point>vertexUp, vector<Point>  vertexDown) {
	vector<Point> res;
	int i = 2;
	while(i<=vertexUp.size())
	{
		int j = i - 1;
		if (i == vertexUp.size()) i = 1;
		float x = (vertexUp[0].x + vertexUp[j].x + vertexUp[i].x) / 3;
		float y = (vertexUp[0].y + vertexUp[j].y + vertexUp[i].y) / 3;
		float z = (vertexUp[0].z + vertexUp[j].z + vertexUp[i].z) / 3;

		res.push_back(Point(x,y,z ));

		x = (vertexDown[j].x + vertexUp[j].x + vertexUp[i].x) / 3;
		y = (vertexDown[j].y + vertexUp[j].y + vertexUp[i].y) / 3;
		z = (vertexDown[j].z + vertexUp[j].z + vertexUp[i].z) / 3;

		res.push_back(Point(x, y, z));

		if (i == 1)break;//完成一周
		i++;
	}

	i = 2;
	while (i <= vertexDown.size())
	{
		int j = i - 1;
		if (i == vertexDown.size()) i = 1;

		float x = (vertexDown[0].x + vertexDown[j].x + vertexDown[i].x) / 3;
		float y = (vertexDown[0].y + vertexDown[j].y + vertexDown[i].y) / 3;
		float z = (vertexDown[0].z + vertexDown[j].z + vertexDown[i].z) / 3;

		res.push_back(Point(x, y, z));

		x = (vertexDown[j].x + vertexUp[i].x + vertexDown[i].x) / 3;
		y = (vertexDown[j].y + vertexUp[i].y + vertexDown[i].y) / 3;
		z = (vertexDown[j].z + vertexUp[i].z + vertexDown[i].z) / 3;

		res.push_back(Point(x, y, z));

		if (i == 1)break;//完成一周

		i++;
	}

	return res;
}

//double angle(Point p1,Point p2) {
//
//	float vectorA[] = { p1.x,p1.y,p1.z };
//	float vectorB[] = { p2.x,p2.y,p2.z };
//
//	double dotProduct = inner_product(vectorA, vectorA + 3, vectorB, 0.0);
//
//	double cos = dotProduct / ( sqrt(pow(p1.x, 2) + pow(p1.y, 2))+pow(p1.z,2)*sqrt(pow(p2.x, 2) + pow(p2.y, 2)) + pow(p2.z, 2));
//
//	double res = acos(cos);
//}