/*****************************************************************************
* | File      	:   LCD_2in_test.c
* | Author      :   Waveshare team
* | Function    :   2inch LCD test demo
* | Info        :
*----------------
* |	This version:   V1.0
* | Date        :   2021-08-20
* | Info        :
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documnetation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to  whom the Software is
# furished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS OR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
******************************************************************************/
#include "LCD_Test.h"
#include "LCD_2in.h"
#include "ImageData.h"
#include "pico/time.h"

bool reserved_addr(uint8_t addr) {
return (addr & 0x78) == 0 || (addr & 0x78) == 0x78;
}

int LCD_2in_test(uint8_t *screen_buff, uint32_t icount)
{

    printf("LCD_2in_test Demo\r\n");
    // if (do_init>0) {
    //     if(DEV_Module_Init()!=0){
    //         return(-1);
    //     }
    // }

    // DEV_SET_PWM(50);
    // /* LCD Init */
    // printf("2inch LCD demo...\r\n");
    // LCD_2IN_Init(HORIZONTAL);
    // LCD_2IN_Clear(WHITE);
    // Paint_SetRotate(ROTATE_90);
    
    //LCD_SetBacklight(1023);
    // UDOUBLE Imagesize = LCD_2IN_HEIGHT*LCD_2IN_WIDTH*2;
    // UWORD *BlackImage;
    // BlackImage = (UWORD *)malloc(Imagesize);

    // if((BlackImage = (UWORD *)malloc(Imagesize)) == NULL) {
    // if (BlackImage == NULL) {
    //     printf("Failed to apply for black memory...\r\n");
    //     return(-2);
    // }
    // /*1.Create a new image cache named IMAGE_RGB and fill it with white*/
    
    // /* GUI */
    printf("drawing...\r\n");
    // /*2.Drawing on the image*/
#if 1
     Paint_Clear(WHITE);
    //  Paint_DrawImage1(gImage_2inch_1,0,0,320,240);
     LCD_2IN_Display((uint8_t *)screen_buff);
    //  DEV_Delay_ms(1000);
     sleep_ms(1000);
#endif

#if 1
     Paint_Clear(BLUE);
    //  Paint_DrawImage1(gImage_2inch_1,0,0,320,240);
     LCD_2IN_Display((uint8_t *)screen_buff);
    //  DEV_Delay_ms(1000);
     sleep_ms(1000);
#endif

#if 1
     Paint_Clear(RED);
    //  Paint_DrawImage1(gImage_2inch_1,0,0,320,240);
     LCD_2IN_Display((uint8_t *)screen_buff);
    //  DEV_Delay_ms(1000);
     sleep_ms(1000);
#endif


#if 1
    Paint_Clear(WHITE);
    Paint_DrawPoint(2,1, BLACK, DOT_PIXEL_1X1,  DOT_FILL_RIGHTUP);//240 240
    Paint_DrawPoint(2,6, BLACK, DOT_PIXEL_2X2,  DOT_FILL_RIGHTUP);
    Paint_DrawPoint(2,11, BLACK, DOT_PIXEL_3X3, DOT_FILL_RIGHTUP);
    Paint_DrawPoint(2,16, BLACK, DOT_PIXEL_4X4, DOT_FILL_RIGHTUP);
    Paint_DrawPoint(2,21, BLACK, DOT_PIXEL_5X5, DOT_FILL_RIGHTUP);
    Paint_DrawLine( 10,  5, 40, 35, MAGENTA, DOT_PIXEL_2X2, LINE_STYLE_SOLID);
    Paint_DrawLine( 10, 35, 40,  5, MAGENTA, DOT_PIXEL_2X2, LINE_STYLE_SOLID);

    Paint_DrawLine( 80,  20, 110, 20, CYAN, DOT_PIXEL_1X1, LINE_STYLE_DOTTED);
    Paint_DrawLine( 95,   5,  95, 35, CYAN, DOT_PIXEL_1X1, LINE_STYLE_DOTTED);

    Paint_DrawRectangle(10, 5, 40, 35, RED, DOT_PIXEL_2X2,DRAW_FILL_EMPTY);
    Paint_DrawRectangle(45, 5, 75, 35, BLUE, DOT_PIXEL_2X2,DRAW_FILL_FULL);

    Paint_DrawCircle(95, 20, 15, GREEN, DOT_PIXEL_1X1, DRAW_FILL_EMPTY);
    Paint_DrawCircle(130, 20, 15, GREEN, DOT_PIXEL_1X1, DRAW_FILL_FULL);


    // Paint_DrawNum (50, 40 ,9.87654321, &Font20,5,  WHITE,  BLACK); // this breaks things 
    //Paint_DrawNum (50, 40 ,9.8, &Font20,1,  WHITE,  BLACK); // this breaks things 

    char buff[256];
    snprintf(buff, sizeof(buff), "%d", icount);
    Paint_DrawString_EN(50, 40, buff, &Font20,  WHITE,  BLACK);

    Paint_DrawString_EN(1, 40, "ABC", &Font20, 0x000f, 0xfff0);
    //Paint_DrawString_CN(1,60, "��ӭʹ��",  &Font24CN, WHITE, BLUE);
    Paint_DrawString_EN(1, 100, "WaveShare", &Font16, RED, WHITE); 

    // /*3.Refresh the picture in RAM to LCD*/
    LCD_2IN_Display((uint8_t *)screen_buff);
    // DEV_Delay_ms(1000);
    sleep_ms(1000);
#endif

#if 0
    Paint_Clear(WHITE);
    
    Paint_DrawString_EN(1, 40, "ABC", &Font20, 0x000f, 0xfff0);
    
    Paint_DrawString_EN(1, 100, "WaveShare", &Font16, RED, WHITE); 

    // /*3.Refresh the picture in RAM to LCD*/
    LCD_2IN_Display((uint8_t *)screen_buff);
    // DEV_Delay_ms(1000);
    sleep_ms(1000);
#endif
	// 4.Test button
//    int key0 = 15; 
//    int key1 = 17; 
//    int key2 = 2; 
//    int key3 = 3; 
   
//    SET_Infrared_PIN(key0);    
//    SET_Infrared_PIN(key1);
//    SET_Infrared_PIN(key2);
//    SET_Infrared_PIN(key3);

    // Paint_Clear(WHITE);
    // LCD_2IN_Display((UBYTE *)BlackImage);

    // LCD_2IN_Clear(WHITE);
    // DEV_Delay_ms(2000);
//    LCD_2IN_Display((uint8_t * )BlackImage);
   
#if 0   
   while(1){
   	Paint_DrawString_EN(70, 100, " Key     Test", &Font20, WHITE, RED);
       if(DEV_Digital_Read(key0 ) == 0){
       
       	Paint_DrawRectangle(288, 208, 308, 228, YELLOW, DOT_PIXEL_1X1,DRAW_FILL_FULL);
		
       }else  {
       	Paint_DrawRectangle(288, 208, 308, 228, RED, DOT_PIXEL_1X1,DRAW_FILL_FULL);
       }
           
       if(DEV_Digital_Read(key1 ) == 0){
       
          Paint_DrawRectangle(12, 208, 32, 228, YELLOW, DOT_PIXEL_1X1,DRAW_FILL_FULL);
			
       }else  {
           
        	Paint_DrawRectangle(12, 208, 32, 228, RED, DOT_PIXEL_1X1,DRAW_FILL_FULL);
       }
       
       if(DEV_Digital_Read(key2) == 0){
       
           
			Paint_DrawRectangle(12, 12, 32, 32, YELLOW, DOT_PIXEL_1X1,DRAW_FILL_FULL);
       }else  {
           
       	Paint_DrawRectangle(12, 12, 32, 32, RED, DOT_PIXEL_1X1,DRAW_FILL_FULL);
       }
       
       if(DEV_Digital_Read(key3 ) == 0){
       
           
            Paint_DrawRectangle(288, 12, 308, 32, YELLOW, DOT_PIXEL_1X1,DRAW_FILL_FULL);

       }else{
           
            Paint_DrawRectangle(288, 12, 308, 32, RED, DOT_PIXEL_1X1,DRAW_FILL_FULL);   
      
       }		
		LCD_2IN_Display((uint8_t * )BlackImage);             
   }
#endif

    /* Module Exit */
    // free(BlackImage);
    // BlackImage = NULL;
    
    
    // DEV_Module_Exit();
    return 0;
}
