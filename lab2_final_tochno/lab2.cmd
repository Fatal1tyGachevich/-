/*****************************************************************
* lab2.cmd
******************************************************************/
MEMORY
{
 IRAM: o = 00000000h l = 00010000h 
 CE0: o = 80000000h l = 01000000h 
}
SECTIONS
{
 .cinit > IRAM
 .text > IRAM
 .stack > IRAM
 .bss > IRAM
 .const > IRAM
 .data > IRAM
 .far > IRAM
 .switch > IRAM
 .sysmem > IRAM
 .tables > IRAM
 .cio > IRAM 
}
