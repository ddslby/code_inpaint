// inpaint.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#include "image.h"
#include "inpainting.h"


/*int _tmain(int argc, _TCHAR* argv[])
{
	inpainting test("test.bmp");
	test.process();
	
	return 0;
}*/

int main()
{
	inpainting test("test.bmp");
	test.process();

	return 0;
}
