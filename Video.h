//
//  YUV4MPEG2 - simple video file format
//  D.P. Mitchell  2018/05/28
//
#pragma once
#include "Image.h"

struct Video {
    HANDLE          hFile;
    unsigned char   *pchBuffer;
    int             m_nHeight;
    int             m_nWidth;
    float           m_fGamma;

    int     NewVideo(char *szFileName, Image &im);
    int     WriteFrame(Image &im);
    int     Close();
    float   SetGamma(float fGamma) { return m_fGamma = fGamma; }

            Video() : m_nHeight(0), m_nWidth(0), pchBuffer(0), m_fGamma(2.0) {}
            ~Video() { Close(); }
};
