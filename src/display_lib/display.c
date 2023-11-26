
#include <stdio.h>
#include "pico/stdlib.h"

#include "fonts.h"
#include "DEV_Config.h"
#include "LCD_2in.h"
#include "GUI_Paint.h"
//#include "LCD_Test.h"
//#include "ImageData.h"
#include "display.h"

static const uint16_t DC   =  8;
static const uint16_t CS   =  9;
static const uint16_t SCK  = 10;
static const uint16_t MOSI = 11;
static const uint16_t RST  = 12;
static const uint16_t BL   = 13;

static const uint16_t height = LCD_2IN_HEIGHT;
static const uint16_t width = LCD_2IN_WIDTH;
static const uint16_t lines = 5;
static const uint16_t lineWidth = LCD_2IN_WIDTH;
static const uint16_t textSize = 3;
static const uint16_t textHeight = 8*3;
static const uint16_t textStart = 200;
static uint16_t ticksSinceScreen = 0;
static uint16_t maxTicksSinceScreen = 0;

static const uint16_t lineHeight = (uint16_t)(height/lines);
static const uint16_t textWidth = width-textStart;
static const uint16_t textHeightOffset = (uint16_t)((lineHeight - textHeight)/2);

#define DGREEN 0x1360
#define DBLUE 0x0070
#define DRED 0x8000

UDOUBLE Imagesize = LCD_2IN_HEIGHT*LCD_2IN_WIDTH*2;
uint8_t *screen_buff;

int display_init()
{
    static char buff[256];
    uint16_t lineStart;

    #ifdef DEBUG_USB
        printf("DISPLAY: {\"init\":0,\"rev_cmd\":%6d,\"rev_count\":%6d,\"speed_cmd\":%6.1f,\"speed_meas\":%6.1f,\"direction_CW\":\"%6s\"},\n",0,0,0.0,0.0,"    CW");
    #endif

    screen_buff = (uint8_t *)malloc(Imagesize);
    if (screen_buff == NULL) {
        printf("Failed to apply for screen buffer memory...\r\n");
        return(-1);
    }

    DEV_Delay_ms(100);
    printf("LCD_2in_test Demo\r\n");
    if(DEV_Module_Init()!=0){
        return -1;
    }
    DEV_SET_PWM(50);
    /* LCD Init */
    printf("2inch LCD demo...\r\n");
    LCD_2IN_Init(HORIZONTAL);
    LCD_2IN_Clear(WHITE);
    Paint_SetRotate(ROTATE_90);
    DEV_Delay_ms(100);

    Paint_NewImage((uint8_t *)screen_buff,LCD_2IN.WIDTH,LCD_2IN.HEIGHT, 90, WHITE);
    Paint_SetScale(65);


    // background
    Paint_Clear(WHITE);

    lineStart = 0*lineHeight+1;
    Paint_DrawRectangle(0,lineStart,lineWidth,lineStart+2*lineHeight-1,DRED, DOT_PIXEL_2X2,DRAW_FILL_FULL);

    lineStart = 2*lineHeight+1;
    Paint_DrawRectangle(0,lineStart,lineWidth,lineStart+2*lineHeight-1,DBLUE, DOT_PIXEL_2X2,DRAW_FILL_FULL);

    lineStart = 4*lineHeight+1;
    Paint_DrawRectangle(0,lineStart,lineWidth,lineStart+lineHeight-1,DGREEN, DOT_PIXEL_2X2,DRAW_FILL_FULL);


    // text labels
    lineStart = 0*lineHeight+1;
    snprintf(buff, sizeof(buff), "Turns Max: ");
    Paint_DrawString_EN(8,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DRED);


    lineStart = 1*lineHeight+1;
    snprintf(buff, sizeof(buff), "Turns: ");
    Paint_DrawString_EN(8,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DRED);


    lineStart = 2*lineHeight+1;
    snprintf(buff, sizeof(buff), "Speed Cmd: ");
    Paint_DrawString_EN(8,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DBLUE);


    lineStart = 3*lineHeight+1;
    snprintf(buff, sizeof(buff), "Speed: ");
    Paint_DrawString_EN(8,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DBLUE);


    lineStart = 4*lineHeight+1;
    snprintf(buff, sizeof(buff), "Direction: ");
    Paint_DrawString_EN(8,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DGREEN);


    LCD_2IN_Display((uint8_t *)screen_buff);

    #ifdef DEBUG_USB
        printf("DISPLAY: {\"init\":1,\"rev_cmd\":%6d,\"rev_count\":%6d,\"speed_cmd\":%6.1f,\"speed_meas\":%6.1f,\"direction_CW\":\"%6s\"},\n",0,0,0.0,0.0,"    CW");
    #endif

}

int display_page(uint32_t rev_cmd, uint32_t rev_count, 
    float speed_cmd, float speed_meas, bool direction_CW)
{

    uint16_t lineStart;
    static char buff[256];

    // background
    lineStart = 0*lineHeight+1;
    Paint_DrawRectangle(textStart,lineStart,textStart+textWidth,lineStart+2*lineHeight-1,DRED, DOT_PIXEL_2X2,DRAW_FILL_FULL);

    lineStart = 2*lineHeight+1;
    Paint_DrawRectangle(textStart,lineStart,textStart+textWidth,lineStart+2*lineHeight-1,DBLUE, DOT_PIXEL_2X2,DRAW_FILL_FULL);

    lineStart = 4*lineHeight+1;
    Paint_DrawRectangle(textStart,lineStart,textStart+textWidth,lineStart+lineHeight-1,DGREEN, DOT_PIXEL_2X2,DRAW_FILL_FULL);

    // text
    lineStart = 0*lineHeight+1;
    snprintf(buff, sizeof(buff), "%6d", rev_cmd);
    Paint_DrawString_EN(textStart,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DRED);

    lineStart = 1*lineHeight+1;
    snprintf(buff, sizeof(buff), "%6d", rev_count);
    Paint_DrawString_EN(textStart,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DRED);

    lineStart = 2*lineHeight+1;
    snprintf(buff, sizeof(buff), "%6.1f", speed_cmd);
    Paint_DrawString_EN(textStart,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DBLUE);

    lineStart = 3*lineHeight+1;
    snprintf(buff, sizeof(buff), "%6.1f", speed_meas);
    Paint_DrawString_EN(textStart,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DBLUE);

    lineStart = 4*lineHeight+1;
    if (direction_CW) {
        snprintf(buff, sizeof(buff), "%s", "    CW");
    } else {
        snprintf(buff, sizeof(buff), "%s", "   CCW");
    }
    Paint_DrawString_EN(textStart,lineStart+textHeightOffset, buff, &Font24,  WHITE,  DGREEN);

    LCD_2IN_Display((uint8_t *)screen_buff);

    #ifdef DEBUG_USB
        printf("DISPLAY: {\"init\":1,\"rev_cmd\":%6d,\"rev_count\":%6d,\"speed_cmd\":%6.1f,\"speed_meas\":%6.1f,\"direction_CW\":\"%6s\"},\n",rev_cmd,rev_count,speed_cmd,speed_meas,buff);
    #endif

    return 0;
}

int display_close()
{
    free(screen_buff);
    return 0;
}
