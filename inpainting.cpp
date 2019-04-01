#include "StdAfx.h"
#include ".\inpainting.h"
#include <memory.h>
#include <math.h>
#define MAX(a, b)  (((a) > (b)) ? (a) : (b)) 
#define MIN(a, b)  (((a) < (b)) ? (a) : (b)) 

//构造函数
inpainting::inpainting(char * name){
	COLORREF temp;//temp是32bit，临时记录像素值(R通道:0-7bit，G通道:8-15bit，B通道:16-23bit)
	m_pImage = new CImage;
	m_bOpen = m_pImage->Load(name);//Load(name);  //载入待修复图像
	if(m_bOpen){//打开成功
		m_pImage->Lock();//Lock / unlock memory.
		m_width = m_pImage->m_BMPwidth;				//宽
		m_height = m_pImage->m_BMPheigth;			//高

		m_mark = new int[m_width * m_height];		//每个像素点的标记
		m_confid = new double[m_width * m_height];	//每个像素点的置信度
		memset(m_confid, 0, m_width * m_height * sizeof(double));//置信度初始化为0
		m_pri = new double[m_width * m_height];		//每个像素点的优先级
		m_gray  = new double[m_width * m_height];	//灰度图像
		m_source = new bool[m_width * m_height];	//放置匹配块空间
		m_color = new COLORREF[m_width * m_height];
		m_r = new double[m_width * m_height];		//开辟R通道存储空间
		m_g = new double[m_width * m_height];		//G
		m_b = new double[m_width * m_height];		//B
	}
	else{
		printf("one or more file is not opened!\n");//打开失败
	}
	//初始化破损区域为原图大小
	m_top = m_height;  
    m_bottom = 0; 
	m_left = m_width;
	m_right = 0;
	
	for(int y = 0; y < m_height; y++){		//第y行(y表示行)
		for(int x = 0; x < m_width; x++){	//第x列(x表示列)
			m_pImage->GetPixel(x, y, temp);			//获取像素值，放在temp中
			m_color[y * m_width + x] = temp;		//彩色图像
			m_r[y * m_width + x] = GetRValue(temp);	//R
			m_g[y * m_width + x] = GetGValue(temp);	//G
			m_b[y * m_width + x] = GetBValue(temp);	//B
		}
	}
}
//析构函数(释放占有的内存空间)
inpainting::~inpainting(void){
	if(m_pImage)delete m_pImage;
	if(m_mark)	delete m_mark;
	if(m_source)delete m_source;
	if(m_color)	delete m_color;
	if(m_r)		delete m_r;
	if(m_g)		delete m_g;
	if(m_b)		delete m_b;
}

//转换成灰度
void inpainting::Convert2Gray(void){
	COLORREF cc;
	double r, g, b;
	for(int y = 0; y < m_height; y++){
		for(int x = 0; x < m_width; x++){
			m_pImage->GetPixel(x, y, cc);
			r = GetRValue(cc);
			g = GetGValue(cc); 
			b = GetBValue(cc);
			//Gray=R*0.114 + G*0.587 + B*0.299
			m_gray[y * m_width + x] = (double)((r * 3735 + g * 19267 + b * 9765) / 32767);
		}
	}
}
//************************************************************************************
bool inpainting::process(void)
{
	char path[200];
	char temp[30];
	Convert2Gray();  //转换成灰度图像
	DrawBoundary();  //第一次绘制边界
	draw_source();   //找到可用于修补的完好块
	memset(m_pri, 0, m_width * m_height * sizeof(double));//优先级初试化为0
	//只有边界才存在优先级，其他位置的优先级都为0
	for(int j = m_top; j <= m_bottom; j++){
		for(int i = m_left; i <= m_right; i++){
			if(m_mark[j * m_width + i] == BOUNDARY){
				m_pri[j * m_width + i] = priority(i,j);//若是边界，则计算其优先级
			}
		}
	}
	int count = 0;
	while(TargetExist())//检测是否还有未修复的点，没有则退出
	{
		count++;
		double max_pri = 0;//最高优先级
		int pri_x, pri_y;//具有最高优先级的点的位置
		for(int j = m_top; j <= m_bottom; j++){
			for(int i = m_left; i <= m_right; i++){
				if((m_mark[j * m_width + i] == BOUNDARY) && (m_pri[j * m_width + i] > max_pri)){//寻找边界上具有最高优先级的点
					pri_x = i;
					pri_y = j;
					max_pri = m_pri[j * m_width + i];
				}
			}
		}
		//printf("pri_x is %d, pri_y is %d, amount is %f\n", pri_x, pri_y, max_pri);
		int patch_x, patch_y;
		PatchTexture(pri_x, pri_y, patch_x, patch_y);//找到与待修复点对应块最相似的块 的中心位置
		update(pri_x, pri_y, patch_x, patch_y, ComputeConfidence(pri_x, pri_y));//修复此块并更新置信度
		UpdateBoundary(pri_x, pri_y);	//更新边界线
		UpdatePri(pri_x, pri_y);		//更新优先级

/**************以下是保存一次修复后的图片***********************/
		strcpy(path, save_path);
		if(count < 10)strcat(path, "000");
		else if(count < 100)strcat(path, "00");
		else if(count < 1000)strcat(path, "0");
		_itoa( count, temp, 10 );
		strcat(path, temp);
		strcat(path, ".bmp");
		m_pImage->Save(path);
		//m_pImage->Save("inpaint_result.bmp");
	}
	return true;
}

void inpainting::DrawBoundary(void){//第一次绘制边界线
	COLORREF color;	
	for(int y = 0; y < m_height; y++){
		for(int x = 0; x < m_width; x++){
			m_pImage->GetPixel(x, y, color);

			if (PAINT_COLOR == color){				//如果像素是绿色，则标记为破损区
				m_mark[y * m_width + x] = TARGET;	//破损区
				m_confid[y * m_width + x] = 0;		//置信度为0
			}
			else{
				m_mark[y * m_width + x] = SOURCE;	//不是绿色，则标记为完好区域
				m_confid[y * m_width + x] = 1;		//置信度为0
			}
		}
	}

	for(int j = 0; j< m_height; j++){//行
	    for(int i = 0; i< m_width; i++){//列
			if(m_mark[j * m_width + i] == TARGET){	//是破损区
				if(i < m_left)	m_left = i;			//调整矩形区域大小
				if(i > m_right)	m_right = i;
				if(j > m_bottom)m_bottom = j;
				if(j < m_top)	m_top = j;
				//若破损点(j,i)的四邻域内存在至少一个完好点或(j,i)位置在图像边缘(最外圈)，则标记为边界
				if(j == m_height - 1 || j == 0 || i == 0 || i == m_width - 1//
					|| m_mark[(j - 1) * m_width + i] == SOURCE	//up
					|| m_mark[j * m_width + i - 1] == SOURCE	//left
					|| m_mark[j * m_width + i + 1] == SOURCE	//right
					|| m_mark[(j + 1) * m_width + i] == SOURCE){//down
						m_mark[j * m_width + i] = BOUNDARY;//边界点
				}
			}
		}
	}
}

double inpainting::priority(int i, int j){
	double confidence, data;
	confidence = ComputeConfidence(i, j);	//置信度
	data = ComputeData(i, j);				//数据项
	return confidence * data;
}

double inpainting::ComputeConfidence(int i, int j) {
	double confidence = 0;
	for(int y = MAX(j - winsize, 0); y <= MIN(j + winsize, m_height - 1); y++){//保证位置合法
		for(int x = MAX(i - winsize, 0); x <= MIN(i + winsize, m_width - 1); x++){
			confidence += m_confid[y * m_width + x];//(j,i)的9*9窗口邻域内完好点的个数
		}
	}
	confidence /= (winsize * 2 + 1) * (winsize * 2 + 1);//置信度 = 完好点个数 / 窗口内像素点个数
	return confidence;

}

double inpainting::ComputeData(int i, int j){
	gradient grad;		//存储最大梯度值的梯度
	gradient temp;		//临时记录梯度
	gradient grad_T;	//等照度线向量
	grad.grad_x = 0;
	grad.grad_y = 0;
	double result;		//数据项计算结果
	double magnitude;	//梯度值
	double max = 0;		//记录梯度最大值
	int x, y;
	//在块中找到最大的梯度值
	for(y = MAX(j - winsize, 0); y <= MIN(j + winsize, m_height - 1); y++){//保证位置合法
		for( x = MAX(i - winsize, 0); x <= MIN(i + winsize, m_width - 1); x++){
			// find the greatest gradient in this patch, this will be the gradient of this pixel(according to "detail paper")
			
			if(m_mark[y * m_width + x] >= 0){//是完好点
			//since I use four neighbors to calculate the gradient, make sure this four neighbors do not touch target region(big jump in gradient)
				//保证四邻域也是完好点
				if(m_mark[y * m_width + (x + 1)] < 0		//right是非完好点
					|| m_mark[y * m_width + (x - 1)] < 0	//left...
					|| m_mark[(y + 1) * m_width + x] < 0	//down...
					|| m_mark[(y - 1) * m_width + x] < 0){	//up...
						continue;
				}
 				temp = GetGradient(x, y);	//计算(x,y)的梯度
				magnitude = temp.grad_x * temp.grad_x + temp.grad_y * temp.grad_y;//梯度值
				if(magnitude > max){		//找到最大的梯度值
					grad.grad_x = temp.grad_x;
					grad.grad_y = temp.grad_y;
					max = magnitude;
				}
			}
		}
	}
	grad_T.grad_x = grad.grad_y;//等照度线(垂直于梯度)
	grad_T.grad_y = -grad.grad_x;

	norm nn = GetNorm(i, j);
	result = (nn.norm_x * grad_T.grad_x) + (nn.norm_y * grad_T.grad_y);//点乘
	result /= 255; //归一化因子255
	result = fabs(result);			
	return result;
}



gradient inpainting::GetGradient(int i, int j)
{
	gradient result;
	result.grad_x = (m_gray[j * m_width + (i + 1)] - m_gray[j * m_width + (i - 1)]) / 2.0;//中心差分
	result.grad_y = (m_gray[(j + 1) * m_width + i] - m_gray[(j - 1) * m_width + i]) / 2.0;

	if(i == 0){//在图像的最外圈边界上(最左边)
		result.grad_x = m_gray[j * m_width + (i + 1)] - m_gray[j * m_width + i];
	}
	if(i == m_width - 1){//最右边
		result.grad_x = m_gray[j * m_width + i] - m_gray[j * m_width + (i - 1)];
	}
	if(j == 0){//最上边
		result.grad_y = m_gray[(j + 1) * m_width + i] - m_gray[j * m_width + i];
	}
	if(j == m_height - 1){//最下边
		result.grad_y = m_gray[j * m_width + i] - m_gray[(j - 1) * m_width + i];
	}
	return result;//返回梯度
}

norm inpainting::GetNorm(int i, int j)
{
	norm result;
	int num = 0;//窗口中边界点的个数
	int neighbor_x[9];//3x3的窗口
	int neighbor_y[9];
	int record[9];
	int count = 0;//窗口总像素个数(只有处于图像的最外圈时，count不是9)
	for(int y = MAX(j - 1, 0); y <= MIN(j + 1, m_height - 1); y++){//保证位置合法
		for(int x = MAX(i - 1, 0); x <= MIN(i + 1, m_width - 1); x++){//扫描3x3的窗口
			count++;
			if(x == i && y == j)continue;	//不包含中心点
			if(m_mark[y * m_width + x] == BOUNDARY){
				num++;			//窗口中边界点的个数
				neighbor_x[num] = x;
				neighbor_y[num] = y;
				record[num] = count;
			}
		}
	}
	if(num == 0 || num == 1){//若不存在2邻域，则随机赋值
		result.norm_x = 0.6;
		result.norm_y = 0.8;
		return result;
	}
	// draw a line between the two neighbors of the boundary pixel, then the norm is the perpendicular to the line
	//划线
	int n_x = neighbor_x[2] - neighbor_x[1];
	int n_y = neighbor_y[2] - neighbor_y[1];
	int temp = n_x;
	n_x = n_y;
	n_y = temp;
	double square = pow(double(n_x * n_x + n_y * n_y), 0.5);//算术平方根号
		
	result.norm_x = n_x / square;
	result.norm_y = n_y / square;
	return result;
}

bool inpainting::draw_source(void)
{
	//在像素周围画一个窗口，若窗口内所有点都是完好的，则这个块可以用于修补
	bool flag;
	for(int j = 0; j < m_height; j++){
		for(int i = 0; i < m_width; i++){
			flag=1;
			if(i < winsize || j < winsize || i >= m_width - winsize || j >= m_height - winsize){
				m_source[j * m_width + i] = false;//无法形成完整的窗口，也就是不能作为
			}
			else{
				for(int y = j - winsize; y <= j + winsize; y++){//扫描窗口内的点是否都是完好点
					for(int x = i - winsize; x <= i + winsize; x++){
						if(m_mark[y * m_width + x] != SOURCE){	//窗口若存在一个非完好点，则结束本窗口的扫描；转而扫描下一个窗口
							m_source[j * m_width + i] = false;
							flag = false;
							break;			
						}		
					}
					if(flag == false)
						break;
				}
			    if(flag != false)
					m_source[j * m_width + i] = true;//若窗口内都是完好点，则标记为可用于修复的块
			}
		}
	}
	return true;
}

bool inpainting::PatchTexture(int x, int y, int &patch_x, int &patch_y)
{
	double temp_r;//R的误差
	double temp_g;//G的误差
	double temp_b;//B的误差
	//根据SSD找到最相似的块(Sum of Squared Differences:误差平方和，计算的是两个块的L2距离)
//	COLORREF color_target, color_source, color_diff;
//	double r0, g0, b0;
//	double r1, g1, b1;

	long min = 99999999;
	long sum;
	int source_x, source_y;
	int target_x, target_y;
	for(int j = 0; j < m_height; j++){
		for(int i = 0; i < m_width; i++){		
			if(m_source[j * m_width + i] == false){//寻找可以用于修复的块
				continue;//没找到，接着找
			}
			sum = 0;
			for(int iter_y = (-1) * winsize; iter_y <= winsize; iter_y++){//窗口内扫描
				for(int iter_x = (-1) * winsize; iter_x <= winsize; iter_x++){
					source_x = i + iter_x;
					source_y = j + iter_y;

					target_x = x + iter_x;
					target_y = y + iter_y;
					
					if(target_x < 0 || target_x >= m_width || target_y < 0 || target_y >= m_height){
						continue;
					}
					if(m_mark[target_y * m_width + target_x] >= 0){//扫描到的点是完好点(也就是说:所有对应的完好点之间作差)
						temp_r = m_r[target_y * m_width + target_x] - m_r[source_y * m_width + source_x];
						temp_g = m_g[target_y * m_width + target_x] - m_g[source_y * m_width + source_x];
						temp_b = m_b[target_y * m_width + target_x] - m_b[source_y * m_width + source_x];
						sum += temp_r * temp_r + temp_g * temp_g + temp_b * temp_b;//平方和
					}
				}
			}
			if(sum < min){//找到更相似的，更新位置和更新SSD值
				min = sum;
				patch_x = i;
				patch_y = j;
			}
		}
	}
	//	printf("patch_x is %d, patch_y is %d\n", patch_x, patch_y);
	return true;
}

bool inpainting::update(int target_x, int target_y, int source_x, int source_y, double confid)
{						//(target_x,target_x):待修复点
						//(source_x,source_x):与待修复点对应块最相似 的块 的中心
//	COLORREF color;
//	double r, g, b;

	int x0, y0, x1, y1;
	for(int iter_y = (-1) * winsize; iter_y <= winsize; iter_y++){//窗口扫描，中心是(source_x,source_y)
		for(int iter_x = (-1) * winsize; iter_x <= winsize; iter_x++){
			x0 = source_x + iter_x;//完好点所在的窗口内的点
			y0 = source_y + iter_y;

			x1 = target_x + iter_x;//待修复点所在窗口内的点
			y1 = target_y + iter_y;
			if(m_mark[y1 * m_width + x1] < 0){//非完好点
				m_pImage->SetPixel(x1, y1, m_color[y0 * m_width + x0]);//将(x1,y1)的颜色拷贝到(x0,y0)处，边界颜色被覆盖，边界往内推移
				m_color[y1 * m_width + x1] = m_color[y0 * m_width + x0];
				m_r[y1 * m_width + x1] = m_r[y0 * m_width + x0];//直接拷贝像素值
				m_g[y1 * m_width + x1] = m_g[y0 * m_width + x0];
				m_b[y1 * m_width + x1] = m_b[y0 * m_width + x0];
				//更新灰度图像(即:修复了)
				m_gray[y1 * m_width + x1] = (double)((m_r[y0 * m_width + x0] * 3735 + m_g[y0 * m_width + x0] * 19267 + m_b[y0 * m_width + x0] * 9765) / 32767);
				m_confid[y1 * m_width + x1] = confid;//更新置信度
			}
		}
	}
	return true;
}

bool inpainting::TargetExist(void){
	for(int j = m_top; j <= m_bottom; j++){
		for(int i = m_left; i <= m_right; i++){
			if(m_mark[j * m_width + i] < 0){
				return true;//存在非完好点
			}
		}
	}
	return false;
}

void inpainting::UpdateBoundary(int i, int j){
	//只更窗口新附近(±2个像素)的区域，即:将窗口的半径扩大成了winsize+2.
	int x, y;
	COLORREF color;//记录颜色

	for(y = MAX(j - winsize - 2, 0); y <= MIN(j + winsize + 2, m_height - 1); y++){//保证位置合法
		for(x = MAX(i - winsize - 2, 0); x <= MIN(i + winsize + 2, m_width - 1); x++){
			m_pImage->GetPixel(x, y, color);		//获取颜色
			if(color == PAINT_COLOR){				
				m_mark[y * m_width + x] = TARGET;	//若是边界(绿色)，则标记为破损点
			}
			else{
				m_mark[y * m_width + x] = SOURCE;	//否则，标记为完好点
			}
		}
	}

	for(y = MAX(j - winsize - 2, 0); y <= MIN(j + winsize + 2, m_height - 1); y++){//保证位置合法
		for(x = MAX(i - winsize - 2, 0); x <= MIN(i + winsize + 2, m_width - 1); x++){
			if(m_mark[y * m_width + x] == TARGET){//若是破损点
				if(y == m_height - 1 || y == 0 || x == 0 || x == m_width - 1//在图像边缘上(最外圈)或4邻域中至少有一个是完好点，则标记为边界
					|| m_mark[(y - 1) * m_width + x] == SOURCE
					|| m_mark[y * m_width + x - 1] == SOURCE
					|| m_mark[y * m_width + x + 1] == SOURCE
					|| m_mark[(y + 1) * m_width + x] == SOURCE){
						m_mark[y * m_width + x] = BOUNDARY;//标记为边界
				}
			}
		}
	}
}

void inpainting::UpdatePri(int i, int j){
	//只更窗口新附近(±3个像素)的区域，即:将窗口的半径扩大成了winsize+3.
	int x, y;
	for(y = MAX(j - winsize - 3, 0); y <= MIN(j + winsize + 3, m_height - 1); y++){
		for(x = MAX(i - winsize - 3, 0); x <= MIN(i + winsize + 3, m_width - 1); x++){
			if(m_mark[y * m_width + x] == BOUNDARY){//若是边界，则更新其优先级
				m_pri[y * m_width + x] = priority(x, y);
			}
		}
	}
}
