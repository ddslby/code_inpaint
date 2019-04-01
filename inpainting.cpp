#include "StdAfx.h"
#include ".\inpainting.h"
#include <memory.h>
#include <math.h>
#define MAX(a, b)  (((a) > (b)) ? (a) : (b)) 
#define MIN(a, b)  (((a) < (b)) ? (a) : (b)) 

//���캯��
inpainting::inpainting(char * name){
	COLORREF temp;//temp��32bit����ʱ��¼����ֵ(Rͨ��:0-7bit��Gͨ��:8-15bit��Bͨ��:16-23bit)
	m_pImage = new CImage;
	m_bOpen = m_pImage->Load(name);//Load(name);  //������޸�ͼ��
	if(m_bOpen){//�򿪳ɹ�
		m_pImage->Lock();//Lock / unlock memory.
		m_width = m_pImage->m_BMPwidth;				//��
		m_height = m_pImage->m_BMPheigth;			//��

		m_mark = new int[m_width * m_height];		//ÿ�����ص�ı��
		m_confid = new double[m_width * m_height];	//ÿ�����ص�����Ŷ�
		memset(m_confid, 0, m_width * m_height * sizeof(double));//���Ŷȳ�ʼ��Ϊ0
		m_pri = new double[m_width * m_height];		//ÿ�����ص�����ȼ�
		m_gray  = new double[m_width * m_height];	//�Ҷ�ͼ��
		m_source = new bool[m_width * m_height];	//����ƥ���ռ�
		m_color = new COLORREF[m_width * m_height];
		m_r = new double[m_width * m_height];		//����Rͨ���洢�ռ�
		m_g = new double[m_width * m_height];		//G
		m_b = new double[m_width * m_height];		//B
	}
	else{
		printf("one or more file is not opened!\n");//��ʧ��
	}
	//��ʼ����������Ϊԭͼ��С
	m_top = m_height;  
    m_bottom = 0; 
	m_left = m_width;
	m_right = 0;
	
	for(int y = 0; y < m_height; y++){		//��y��(y��ʾ��)
		for(int x = 0; x < m_width; x++){	//��x��(x��ʾ��)
			m_pImage->GetPixel(x, y, temp);			//��ȡ����ֵ������temp��
			m_color[y * m_width + x] = temp;		//��ɫͼ��
			m_r[y * m_width + x] = GetRValue(temp);	//R
			m_g[y * m_width + x] = GetGValue(temp);	//G
			m_b[y * m_width + x] = GetBValue(temp);	//B
		}
	}
}
//��������(�ͷ�ռ�е��ڴ�ռ�)
inpainting::~inpainting(void){
	if(m_pImage)delete m_pImage;
	if(m_mark)	delete m_mark;
	if(m_source)delete m_source;
	if(m_color)	delete m_color;
	if(m_r)		delete m_r;
	if(m_g)		delete m_g;
	if(m_b)		delete m_b;
}

//ת���ɻҶ�
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
	Convert2Gray();  //ת���ɻҶ�ͼ��
	DrawBoundary();  //��һ�λ��Ʊ߽�
	draw_source();   //�ҵ��������޲�����ÿ�
	memset(m_pri, 0, m_width * m_height * sizeof(double));//���ȼ����Ի�Ϊ0
	//ֻ�б߽�Ŵ������ȼ�������λ�õ����ȼ���Ϊ0
	for(int j = m_top; j <= m_bottom; j++){
		for(int i = m_left; i <= m_right; i++){
			if(m_mark[j * m_width + i] == BOUNDARY){
				m_pri[j * m_width + i] = priority(i,j);//���Ǳ߽磬����������ȼ�
			}
		}
	}
	int count = 0;
	while(TargetExist())//����Ƿ���δ�޸��ĵ㣬û�����˳�
	{
		count++;
		double max_pri = 0;//������ȼ�
		int pri_x, pri_y;//����������ȼ��ĵ��λ��
		for(int j = m_top; j <= m_bottom; j++){
			for(int i = m_left; i <= m_right; i++){
				if((m_mark[j * m_width + i] == BOUNDARY) && (m_pri[j * m_width + i] > max_pri)){//Ѱ�ұ߽��Ͼ���������ȼ��ĵ�
					pri_x = i;
					pri_y = j;
					max_pri = m_pri[j * m_width + i];
				}
			}
		}
		//printf("pri_x is %d, pri_y is %d, amount is %f\n", pri_x, pri_y, max_pri);
		int patch_x, patch_y;
		PatchTexture(pri_x, pri_y, patch_x, patch_y);//�ҵ�����޸����Ӧ�������ƵĿ� ������λ��
		update(pri_x, pri_y, patch_x, patch_y, ComputeConfidence(pri_x, pri_y));//�޸��˿鲢�������Ŷ�
		UpdateBoundary(pri_x, pri_y);	//���±߽���
		UpdatePri(pri_x, pri_y);		//�������ȼ�

/**************�����Ǳ���һ���޸����ͼƬ***********************/
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

void inpainting::DrawBoundary(void){//��һ�λ��Ʊ߽���
	COLORREF color;	
	for(int y = 0; y < m_height; y++){
		for(int x = 0; x < m_width; x++){
			m_pImage->GetPixel(x, y, color);

			if (PAINT_COLOR == color){				//�����������ɫ������Ϊ������
				m_mark[y * m_width + x] = TARGET;	//������
				m_confid[y * m_width + x] = 0;		//���Ŷ�Ϊ0
			}
			else{
				m_mark[y * m_width + x] = SOURCE;	//������ɫ������Ϊ�������
				m_confid[y * m_width + x] = 1;		//���Ŷ�Ϊ0
			}
		}
	}

	for(int j = 0; j< m_height; j++){//��
	    for(int i = 0; i< m_width; i++){//��
			if(m_mark[j * m_width + i] == TARGET){	//��������
				if(i < m_left)	m_left = i;			//�������������С
				if(i > m_right)	m_right = i;
				if(j > m_bottom)m_bottom = j;
				if(j < m_top)	m_top = j;
				//�������(j,i)���������ڴ�������һ����õ��(j,i)λ����ͼ���Ե(����Ȧ)������Ϊ�߽�
				if(j == m_height - 1 || j == 0 || i == 0 || i == m_width - 1//
					|| m_mark[(j - 1) * m_width + i] == SOURCE	//up
					|| m_mark[j * m_width + i - 1] == SOURCE	//left
					|| m_mark[j * m_width + i + 1] == SOURCE	//right
					|| m_mark[(j + 1) * m_width + i] == SOURCE){//down
						m_mark[j * m_width + i] = BOUNDARY;//�߽��
				}
			}
		}
	}
}

double inpainting::priority(int i, int j){
	double confidence, data;
	confidence = ComputeConfidence(i, j);	//���Ŷ�
	data = ComputeData(i, j);				//������
	return confidence * data;
}

double inpainting::ComputeConfidence(int i, int j) {
	double confidence = 0;
	for(int y = MAX(j - winsize, 0); y <= MIN(j + winsize, m_height - 1); y++){//��֤λ�úϷ�
		for(int x = MAX(i - winsize, 0); x <= MIN(i + winsize, m_width - 1); x++){
			confidence += m_confid[y * m_width + x];//(j,i)��9*9������������õ�ĸ���
		}
	}
	confidence /= (winsize * 2 + 1) * (winsize * 2 + 1);//���Ŷ� = ��õ���� / ���������ص����
	return confidence;

}

double inpainting::ComputeData(int i, int j){
	gradient grad;		//�洢����ݶ�ֵ���ݶ�
	gradient temp;		//��ʱ��¼�ݶ�
	gradient grad_T;	//���ն�������
	grad.grad_x = 0;
	grad.grad_y = 0;
	double result;		//�����������
	double magnitude;	//�ݶ�ֵ
	double max = 0;		//��¼�ݶ����ֵ
	int x, y;
	//�ڿ����ҵ������ݶ�ֵ
	for(y = MAX(j - winsize, 0); y <= MIN(j + winsize, m_height - 1); y++){//��֤λ�úϷ�
		for( x = MAX(i - winsize, 0); x <= MIN(i + winsize, m_width - 1); x++){
			// find the greatest gradient in this patch, this will be the gradient of this pixel(according to "detail paper")
			
			if(m_mark[y * m_width + x] >= 0){//����õ�
			//since I use four neighbors to calculate the gradient, make sure this four neighbors do not touch target region(big jump in gradient)
				//��֤������Ҳ����õ�
				if(m_mark[y * m_width + (x + 1)] < 0		//right�Ƿ���õ�
					|| m_mark[y * m_width + (x - 1)] < 0	//left...
					|| m_mark[(y + 1) * m_width + x] < 0	//down...
					|| m_mark[(y - 1) * m_width + x] < 0){	//up...
						continue;
				}
 				temp = GetGradient(x, y);	//����(x,y)���ݶ�
				magnitude = temp.grad_x * temp.grad_x + temp.grad_y * temp.grad_y;//�ݶ�ֵ
				if(magnitude > max){		//�ҵ������ݶ�ֵ
					grad.grad_x = temp.grad_x;
					grad.grad_y = temp.grad_y;
					max = magnitude;
				}
			}
		}
	}
	grad_T.grad_x = grad.grad_y;//���ն���(��ֱ���ݶ�)
	grad_T.grad_y = -grad.grad_x;

	norm nn = GetNorm(i, j);
	result = (nn.norm_x * grad_T.grad_x) + (nn.norm_y * grad_T.grad_y);//���
	result /= 255; //��һ������255
	result = fabs(result);			
	return result;
}



gradient inpainting::GetGradient(int i, int j)
{
	gradient result;
	result.grad_x = (m_gray[j * m_width + (i + 1)] - m_gray[j * m_width + (i - 1)]) / 2.0;//���Ĳ��
	result.grad_y = (m_gray[(j + 1) * m_width + i] - m_gray[(j - 1) * m_width + i]) / 2.0;

	if(i == 0){//��ͼ�������Ȧ�߽���(�����)
		result.grad_x = m_gray[j * m_width + (i + 1)] - m_gray[j * m_width + i];
	}
	if(i == m_width - 1){//���ұ�
		result.grad_x = m_gray[j * m_width + i] - m_gray[j * m_width + (i - 1)];
	}
	if(j == 0){//���ϱ�
		result.grad_y = m_gray[(j + 1) * m_width + i] - m_gray[j * m_width + i];
	}
	if(j == m_height - 1){//���±�
		result.grad_y = m_gray[j * m_width + i] - m_gray[(j - 1) * m_width + i];
	}
	return result;//�����ݶ�
}

norm inpainting::GetNorm(int i, int j)
{
	norm result;
	int num = 0;//�����б߽��ĸ���
	int neighbor_x[9];//3x3�Ĵ���
	int neighbor_y[9];
	int record[9];
	int count = 0;//���������ظ���(ֻ�д���ͼ�������Ȧʱ��count����9)
	for(int y = MAX(j - 1, 0); y <= MIN(j + 1, m_height - 1); y++){//��֤λ�úϷ�
		for(int x = MAX(i - 1, 0); x <= MIN(i + 1, m_width - 1); x++){//ɨ��3x3�Ĵ���
			count++;
			if(x == i && y == j)continue;	//���������ĵ�
			if(m_mark[y * m_width + x] == BOUNDARY){
				num++;			//�����б߽��ĸ���
				neighbor_x[num] = x;
				neighbor_y[num] = y;
				record[num] = count;
			}
		}
	}
	if(num == 0 || num == 1){//��������2�����������ֵ
		result.norm_x = 0.6;
		result.norm_y = 0.8;
		return result;
	}
	// draw a line between the two neighbors of the boundary pixel, then the norm is the perpendicular to the line
	//����
	int n_x = neighbor_x[2] - neighbor_x[1];
	int n_y = neighbor_y[2] - neighbor_y[1];
	int temp = n_x;
	n_x = n_y;
	n_y = temp;
	double square = pow(double(n_x * n_x + n_y * n_y), 0.5);//����ƽ������
		
	result.norm_x = n_x / square;
	result.norm_y = n_y / square;
	return result;
}

bool inpainting::draw_source(void)
{
	//��������Χ��һ�����ڣ������������е㶼����õģ����������������޲�
	bool flag;
	for(int j = 0; j < m_height; j++){
		for(int i = 0; i < m_width; i++){
			flag=1;
			if(i < winsize || j < winsize || i >= m_width - winsize || j >= m_height - winsize){
				m_source[j * m_width + i] = false;//�޷��γ������Ĵ��ڣ�Ҳ���ǲ�����Ϊ
			}
			else{
				for(int y = j - winsize; y <= j + winsize; y++){//ɨ�贰���ڵĵ��Ƿ�����õ�
					for(int x = i - winsize; x <= i + winsize; x++){
						if(m_mark[y * m_width + x] != SOURCE){	//����������һ������õ㣬����������ڵ�ɨ�裻ת��ɨ����һ������
							m_source[j * m_width + i] = false;
							flag = false;
							break;			
						}		
					}
					if(flag == false)
						break;
				}
			    if(flag != false)
					m_source[j * m_width + i] = true;//�������ڶ�����õ㣬����Ϊ�������޸��Ŀ�
			}
		}
	}
	return true;
}

bool inpainting::PatchTexture(int x, int y, int &patch_x, int &patch_y)
{
	double temp_r;//R�����
	double temp_g;//G�����
	double temp_b;//B�����
	//����SSD�ҵ������ƵĿ�(Sum of Squared Differences:���ƽ���ͣ���������������L2����)
//	COLORREF color_target, color_source, color_diff;
//	double r0, g0, b0;
//	double r1, g1, b1;

	long min = 99999999;
	long sum;
	int source_x, source_y;
	int target_x, target_y;
	for(int j = 0; j < m_height; j++){
		for(int i = 0; i < m_width; i++){		
			if(m_source[j * m_width + i] == false){//Ѱ�ҿ��������޸��Ŀ�
				continue;//û�ҵ���������
			}
			sum = 0;
			for(int iter_y = (-1) * winsize; iter_y <= winsize; iter_y++){//������ɨ��
				for(int iter_x = (-1) * winsize; iter_x <= winsize; iter_x++){
					source_x = i + iter_x;
					source_y = j + iter_y;

					target_x = x + iter_x;
					target_y = y + iter_y;
					
					if(target_x < 0 || target_x >= m_width || target_y < 0 || target_y >= m_height){
						continue;
					}
					if(m_mark[target_y * m_width + target_x] >= 0){//ɨ�赽�ĵ�����õ�(Ҳ����˵:���ж�Ӧ����õ�֮������)
						temp_r = m_r[target_y * m_width + target_x] - m_r[source_y * m_width + source_x];
						temp_g = m_g[target_y * m_width + target_x] - m_g[source_y * m_width + source_x];
						temp_b = m_b[target_y * m_width + target_x] - m_b[source_y * m_width + source_x];
						sum += temp_r * temp_r + temp_g * temp_g + temp_b * temp_b;//ƽ����
					}
				}
			}
			if(sum < min){//�ҵ������Ƶģ�����λ�ú͸���SSDֵ
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
{						//(target_x,target_x):���޸���
						//(source_x,source_x):����޸����Ӧ�������� �Ŀ� ������
//	COLORREF color;
//	double r, g, b;

	int x0, y0, x1, y1;
	for(int iter_y = (-1) * winsize; iter_y <= winsize; iter_y++){//����ɨ�裬������(source_x,source_y)
		for(int iter_x = (-1) * winsize; iter_x <= winsize; iter_x++){
			x0 = source_x + iter_x;//��õ����ڵĴ����ڵĵ�
			y0 = source_y + iter_y;

			x1 = target_x + iter_x;//���޸������ڴ����ڵĵ�
			y1 = target_y + iter_y;
			if(m_mark[y1 * m_width + x1] < 0){//����õ�
				m_pImage->SetPixel(x1, y1, m_color[y0 * m_width + x0]);//��(x1,y1)����ɫ������(x0,y0)�����߽���ɫ�����ǣ��߽���������
				m_color[y1 * m_width + x1] = m_color[y0 * m_width + x0];
				m_r[y1 * m_width + x1] = m_r[y0 * m_width + x0];//ֱ�ӿ�������ֵ
				m_g[y1 * m_width + x1] = m_g[y0 * m_width + x0];
				m_b[y1 * m_width + x1] = m_b[y0 * m_width + x0];
				//���»Ҷ�ͼ��(��:�޸���)
				m_gray[y1 * m_width + x1] = (double)((m_r[y0 * m_width + x0] * 3735 + m_g[y0 * m_width + x0] * 19267 + m_b[y0 * m_width + x0] * 9765) / 32767);
				m_confid[y1 * m_width + x1] = confid;//�������Ŷ�
			}
		}
	}
	return true;
}

bool inpainting::TargetExist(void){
	for(int j = m_top; j <= m_bottom; j++){
		for(int i = m_left; i <= m_right; i++){
			if(m_mark[j * m_width + i] < 0){
				return true;//���ڷ���õ�
			}
		}
	}
	return false;
}

void inpainting::UpdateBoundary(int i, int j){
	//ֻ�������¸���(��2������)�����򣬼�:�����ڵİ뾶�������winsize+2.
	int x, y;
	COLORREF color;//��¼��ɫ

	for(y = MAX(j - winsize - 2, 0); y <= MIN(j + winsize + 2, m_height - 1); y++){//��֤λ�úϷ�
		for(x = MAX(i - winsize - 2, 0); x <= MIN(i + winsize + 2, m_width - 1); x++){
			m_pImage->GetPixel(x, y, color);		//��ȡ��ɫ
			if(color == PAINT_COLOR){				
				m_mark[y * m_width + x] = TARGET;	//���Ǳ߽�(��ɫ)������Ϊ�����
			}
			else{
				m_mark[y * m_width + x] = SOURCE;	//���򣬱��Ϊ��õ�
			}
		}
	}

	for(y = MAX(j - winsize - 2, 0); y <= MIN(j + winsize + 2, m_height - 1); y++){//��֤λ�úϷ�
		for(x = MAX(i - winsize - 2, 0); x <= MIN(i + winsize + 2, m_width - 1); x++){
			if(m_mark[y * m_width + x] == TARGET){//���������
				if(y == m_height - 1 || y == 0 || x == 0 || x == m_width - 1//��ͼ���Ե��(����Ȧ)��4������������һ������õ㣬����Ϊ�߽�
					|| m_mark[(y - 1) * m_width + x] == SOURCE
					|| m_mark[y * m_width + x - 1] == SOURCE
					|| m_mark[y * m_width + x + 1] == SOURCE
					|| m_mark[(y + 1) * m_width + x] == SOURCE){
						m_mark[y * m_width + x] = BOUNDARY;//���Ϊ�߽�
				}
			}
		}
	}
}

void inpainting::UpdatePri(int i, int j){
	//ֻ�������¸���(��3������)�����򣬼�:�����ڵİ뾶�������winsize+3.
	int x, y;
	for(y = MAX(j - winsize - 3, 0); y <= MIN(j + winsize + 3, m_height - 1); y++){
		for(x = MAX(i - winsize - 3, 0); x <= MIN(i + winsize + 3, m_width - 1); x++){
			if(m_mark[y * m_width + x] == BOUNDARY){//���Ǳ߽磬����������ȼ�
				m_pri[y * m_width + x] = priority(x, y);
			}
		}
	}
}
