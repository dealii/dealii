/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

// common Windows parts
#include "winvideo.h"

#include <dxsdkver.h>
#if _DXSDK_PRODUCT_MAJOR >= 9
// new implementation based on Direct2D
#include "d2dvideo.cpp"
#else // _DXSDK_PRODUCT_MAJOR >= 9

// and another headers
#include <cassert>
#include <stdio.h>
#include <ddraw.h>

#pragma comment(lib, "ddraw.lib")
#pragma comment(lib, "dxguid.lib")

LPDIRECTDRAW7               g_pDD = NULL;        // DirectDraw object
LPDIRECTDRAWSURFACE7        g_pDDSPrimary = NULL;// DirectDraw primary surface
LPDIRECTDRAWSURFACE7        g_pDDSBack = NULL;   // DirectDraw back surface
LPDIRECTDRAWSURFACE7        g_pDDSOverlay = NULL;// DirectDraw overlay surface
LPDIRECTDRAWCLIPPER         g_pClipper = NULL;   // DirectDraw clipping struct
DDOVERLAYFX                 g_OverlayFX;         // DirectDraw overlay effects struct
DDCAPS                      g_DDCaps;            // DirectDraw hardware capabilities struct
DWORD                       g_OverlayFlags = 0;  // DirectDraw overlay flags variable
DWORD                       g_dwXRatio,
                            g_dwYRatio;          // The ratios between the src and dst rects
RECT                        g_rcSrc = {0, 0, 0, 0},
                            g_rcDst = {0, 0, 0, 0};
HANDLE                      g_hVSync;

// check for new DX SDK (8 & 9)
#ifdef DDSCAPS_PRIMARYSURFACELEFT
#include <dxerr8.h>
#pragma comment(lib, "dxerr8.lib")
#else
// old SDK (7)
#include <d3dx.h>
#pragma comment(lib, "d3dx.lib")
#endif

//! Create a dialog box and tell the user what went wrong
bool DisplayError(LPSTR lpstrErr, HRESULT hres)
{
    static bool InError = false;
    int retval = 0;
    if (!InError)
    {
        InError = true;
#ifdef DDSCAPS_PRIMARYSURFACELEFT
        const char *message = hres?DXGetErrorString8A(hres):0;
#else
        char message[256]; if(hres) D3DXGetErrorString(hres, 256, message);
#endif
        retval = MessageBoxA(g_hAppWnd, lpstrErr, hres?message:"Error!", MB_OK|MB_ICONERROR);
        InError = false;
    }
    return false;
}

//! Releases the overlay surface
void DestroyOverlay()
{
    if (g_pClipper)
        g_pClipper->Release();
    if (g_pDDSOverlay) {
        g_pImg = 0; LPDIRECTDRAWSURFACE7 pDDSOverlay(g_pDDSOverlay);
        g_pDDSOverlay = NULL;
        YIELD_TO_THREAD();
        pDDSOverlay->Release(); // be sure nobody uses old value
    }
}

//! Releases the primary surface
void DestroyPrimary()
{
    if (g_pDDSPrimary)
    {
        g_pDDSPrimary->Release();
        g_pDDSPrimary = NULL;
    }
}

//! Releases core DirectDraw objects
void DestroyDDraw()
{
    DestroyPrimary();
    // Release the DDraw object
    if (g_pDD) {
        LPDIRECTDRAW7 pDD(g_pDD); // be sure nobody uses old value
        g_pDD = NULL; Sleep(1); pDD->Release();
    }
}

//! Checks and corrects all boundries for alignment and stretching
void CheckBoundries(void)
{
    // Make sure the coordinates fulfill the stretching requirements.  Often
    // the hardware will require a certain ammount of stretching to do
    // overlays. This stretch factor is held in dwMinOverlayStretch as the
    // stretch factor multiplied by 1000 (to keep an accuracy of 3 decimal places).
    if ((g_DDCaps.dwCaps & DDCAPS_OVERLAYSTRETCH) && (g_DDCaps.dwMinOverlayStretch)
        && (g_dwXRatio < g_DDCaps.dwMinOverlayStretch))
    {
        g_rcDst.right = 2 * GetSystemMetrics(SM_CXSIZEFRAME) + g_rcDst.left + (g_sizex
                                 * (g_DDCaps.dwMinOverlayStretch + 1)) / 1000;
        SetWindowTextA(g_hAppWnd, "Window is too small!");
    }
    else if ((g_DDCaps.dwCaps & DDCAPS_OVERLAYSTRETCH) && (g_DDCaps.dwMaxOverlayStretch)
        && (g_dwXRatio > g_DDCaps.dwMaxOverlayStretch))
    {
        g_rcDst.right = 2 * GetSystemMetrics(SM_CXSIZEFRAME) + g_rcDst.left + (g_sizey
                               * (g_DDCaps.dwMaxOverlayStretch + 999)) / 1000;
        SetWindowTextA(g_hAppWnd, "Window is too large!");
    }
    else if(!g_video->calc_fps) SetWindowText(g_hAppWnd, g_video->title);

    // Recalculate the ratio's for the upcoming calculations
    g_dwXRatio = (g_rcDst.right - g_rcDst.left) * 1000 / (g_rcSrc.right - g_rcSrc.left);
    g_dwYRatio = (g_rcDst.bottom - g_rcDst.top) * 1000 / (g_rcSrc.bottom - g_rcSrc.top);

    // Check to make sure we're within the screen's boundries, if not then fix
    // the problem by adjusting the source rectangle which we draw from.
    if (g_rcDst.left < 0)
    {
        g_rcSrc.left = -g_rcDst.left * 1000 / g_dwXRatio;
        g_rcDst.left = 0;
    }
    if (g_rcDst.right > GetSystemMetrics(SM_CXSCREEN))
    {
        g_rcSrc.right = g_sizex - ((g_rcDst.right - GetSystemMetrics(SM_CXSCREEN)) * 1000 / g_dwXRatio);
        g_rcDst.right = GetSystemMetrics(SM_CXSCREEN);
    }
    if (g_rcDst.bottom > GetSystemMetrics(SM_CYSCREEN))
    {
        g_rcSrc.bottom = g_sizey - ((g_rcDst.bottom - GetSystemMetrics(SM_CYSCREEN)) * 1000 / g_dwYRatio);
        g_rcDst.bottom = GetSystemMetrics(SM_CYSCREEN);
    }
    // I don't know how useful this is... but just in case someone can do it - here's the check.
    if (g_rcDst.top < 0)
    {
        g_rcSrc.top = -g_rcDst.top * 1000 / g_dwYRatio;
        g_rcDst.top = 0;
    }

    // Make sure the coordinates fulfill the alignment requirements
    // these expressions (x & -y) just do alignment by dropping low order bits...
    // so to round up, we add first, then truncate.
    if ((g_DDCaps.dwCaps & DDCAPS_ALIGNBOUNDARYSRC) && g_DDCaps.dwAlignBoundarySrc)
        g_rcSrc.left = (g_rcSrc.left + g_DDCaps.dwAlignBoundarySrc / 2) & -(signed)
            (g_DDCaps.dwAlignBoundarySrc);
    if ((g_DDCaps.dwCaps & DDCAPS_ALIGNSIZESRC) && g_DDCaps.dwAlignSizeSrc)
        g_rcSrc.right = g_rcSrc.left + (g_rcSrc.right - g_rcSrc.left + g_DDCaps.dwAlignSizeSrc
                                   / 2) & -(signed) (g_DDCaps.dwAlignSizeSrc);
    if ((g_DDCaps.dwCaps & DDCAPS_ALIGNBOUNDARYDEST) && g_DDCaps.dwAlignBoundaryDest)
        g_rcDst.left = (g_rcDst.left + g_DDCaps.dwAlignBoundaryDest / 2) & -(signed)
            (g_DDCaps.dwAlignBoundaryDest);
    if ((g_DDCaps.dwCaps & DDCAPS_ALIGNSIZEDEST) && g_DDCaps.dwAlignSizeDest)
        g_rcDst.right = g_rcDst.left + (g_rcDst.right - g_rcDst.left) & -(signed) (g_DDCaps.dwAlignSizeDest);
}

//! Get translated by system color value
DWORD DDColorMatch(IDirectDrawSurface7 * pdds, COLORREF rgb)
{
    COLORREF       rgbT;
    HDC            hdc;
    DWORD          dw = CLR_INVALID;
    DDSURFACEDESC2 ddsd;
    HRESULT        hres;

    //  Use GDI SetPixel to color match for us
    if (rgb != CLR_INVALID && pdds->GetDC(&hdc) == DD_OK) {
        rgbT = GetPixel(hdc, 0, 0);     // Save current pixel value
        SetPixel(hdc, 0, 0, rgb);       // Set our value
        pdds->ReleaseDC(hdc);
    }
    // Now lock the surface so we can read back the converted color
    ddsd.dwSize = sizeof(ddsd);
    while ((hres = pdds->Lock(NULL, &ddsd, 0, NULL)) == DDERR_WASSTILLDRAWING)
        YIELD_TO_THREAD();
    if (hres == DD_OK) {
        dw = *(DWORD *) ddsd.lpSurface;                 // Get DWORD
        if (ddsd.ddpfPixelFormat.dwRGBBitCount < 32)
            dw &= (1 << ddsd.ddpfPixelFormat.dwRGBBitCount) - 1;  // Mask it to bpp
        pdds->Unlock(NULL);
    }
    else return DisplayError("Can't lock primary surface", hres);
    //  Now put the color that was there back.
    if (rgb != CLR_INVALID && pdds->GetDC(&hdc) == DD_OK) {
        SetPixel(hdc, 0, 0, rgbT);
        pdds->ReleaseDC(hdc);
    }
    return dw;
}

//! Load the bitmap and copy it to the overlay surface
bool DrawOverlay()
{
    HRESULT        hRet;       // This is where we put return values from DirectDraw.
    DDSURFACEDESC2 surfDesc;
    // Setup structure
    memset(&surfDesc, 0, sizeof(surfDesc)); surfDesc.dwSize = sizeof(surfDesc);

    hRet = g_pDDSOverlay->Lock(NULL, &surfDesc, DDLOCK_SURFACEMEMORYPTR | DDLOCK_NOSYSLOCK | DDLOCK_WRITEONLY, NULL);
    if (hRet != DD_OK ||  surfDesc.lpSurface == NULL)
        return DisplayError("Can't lock overlay surface", hRet);
    else {
        g_pImg = (unsigned int *)surfDesc.lpSurface;
        //g_pDDSOverlay->Unlock(NULL); is not needed?
    }
    // Setup effects structure
    memset(&g_OverlayFX, 0, sizeof(g_OverlayFX)); g_OverlayFX.dwSize = sizeof(g_OverlayFX);
    // Setup overlay flags.
    g_OverlayFlags = DDOVER_SHOW;
    // Check for destination color keying capability
    if ((g_DDCaps.dwCKeyCaps & DDCKEYCAPS_DESTOVERLAY) && ((g_DDCaps.dwCaps & DDCAPS_OVERLAYCANTCLIP) || (g_DDCaps.dwCKeyCaps & DDCKEYCAPS_NOCOSTOVERLAY) ))
    {
        // If so, we'll use it to clip the bitmap when other windows go on top
        // of us. Just for the record - this color range for color keying (the
        // high/low values) are not heavily supported right now, so for almost
        // all cards, just use the same color for both.
        g_OverlayFX.dckDestColorkey.dwColorSpaceLowValue =
        g_OverlayFX.dckDestColorkey.dwColorSpaceHighValue = DDColorMatch(g_pDDSPrimary, RGBKEY);
        g_OverlayFlags |= DDOVER_DDFX | DDOVER_KEYDESTOVERRIDE;
    } else {
        // If not, we'll setup a clipper for the window.  This will fix the
        // problem on a few video cards - but the ones that don't shouldn't care.
        hRet = g_pDD->CreateClipper(0, &g_pClipper, NULL);
        if (hRet != DD_OK)
            return DisplayError("Can't create clipper", hRet);
        hRet = g_pClipper->SetHWnd(0, g_hAppWnd);
        if (hRet != DD_OK)
            return DisplayError("Can't attach clipper", hRet);
        hRet = g_pDDSPrimary->SetClipper(g_pClipper);
        if (hRet != DD_OK)
            return DisplayError("Can't set clipper", hRet);
    }
    return true;
}

//! Init the primary surface
bool DDPrimaryInit()
{
    HRESULT        hRet;
    DDSURFACEDESC2 ddsd;  // A surface description structure

    // Create the primary surface.  The primary surface is the full screen -
    // since we're a windowed app - we'll just write to the portion of the
    // screen within our window.
    memset(&ddsd, 0, sizeof(ddsd)); // Set all fields of struct to 0 and set .dwSize to
    ddsd.dwSize = sizeof(ddsd);     // Sizeof the variable - these two steps required for most DDraw structs
    ddsd.dwFlags = DDSD_CAPS;       // Set flags for variables we're using...
    ddsd.ddsCaps.dwCaps = DDSCAPS_PRIMARYSURFACE;  // Set the variables we said we would in dwFlags
    hRet = g_pDD->CreateSurface(&ddsd, &g_pDDSPrimary, NULL);
    if (hRet != DD_OK)
        return DisplayError("Can't create primary surface", hRet);
    return true;
}

//! Init DirectDraw Stuff
bool DDInit()
{
    HRESULT hRet;
    g_rcSrc.right = g_sizex;
    g_rcSrc.bottom = g_sizey;

    hRet = DirectDrawCreateEx(NULL, (VOID**)&g_pDD, IID_IDirectDraw7, NULL);
    if (hRet != DD_OK)
        return DisplayError("Can't create DirectDraw7 instance", hRet);

    // Set cooperation level with other windows to be normal (ie. not full screen)
    // You MUST set the cooperation level to be SOMETHING, for windowed apps use
    // DDSCL_NORMAL, for full screen use: DDSCL_EXCLUSIVE | DDSCL_FULLSCREEN.
    hRet = g_pDD->SetCooperativeLevel(g_hAppWnd, DDSCL_NORMAL);
    if (hRet != DD_OK)
        return DisplayError("Can't set cooperative level", hRet);
    return DDPrimaryInit();
}

//! Setup the overlay object
bool DDOverlayInit()
{
    // Get hardware's CAPabilitieS
    memset(&g_DDCaps, 0, sizeof(g_DDCaps));
    g_DDCaps.dwSize = sizeof(g_DDCaps);
    if (g_pDD->GetCaps(&g_DDCaps, 0))
        return DisplayError("Can't get capabilities");

    // Make sure it supports overlays
    if (!(g_DDCaps.dwCaps & DDCAPS_OVERLAY))
        return DisplayError("Hardware doesn't support overlays");

    //DO NOT Make sure it supports stretching (scaling)
    //if (!(g_DDCaps.dwCaps & DDCAPS_OVERLAYSTRETCH)) return false;

    DDSURFACEDESC2              ddsd;  // DirectDraw surface descriptor
    HRESULT                     hRet;  // I'm not even going to try...
    // The pixel formats that we want the surface to be in
    DDPIXELFORMAT               ddpfOverlayFormats[] = {
        {sizeof(DDPIXELFORMAT), DDPF_RGB, 0, 32, 0xFF0000, 0x0FF00, 0x0000FF, 0}, // 32-bit RGB
        {sizeof(DDPIXELFORMAT), DDPF_RGB, 0, 16, 0x007C00, 0x003e0, 0x00001F, 0}, // 16-bit RGB 5:5:5
        {sizeof(DDPIXELFORMAT), DDPF_RGB, 0, 16, 0x00F800, 0x007e0, 0x00001F, 0}, // 16-bit RGB 5:6:5
        {sizeof(DDPIXELFORMAT), DDPF_FOURCC, mmioFOURCC('U','Y','V','Y'), 16, 0, 0, 0, 0}, // UYVY
        {sizeof(DDPIXELFORMAT), DDPF_FOURCC, mmioFOURCC('Y','4','2','2'), 16, 0, 0, 0, 0}, // the same as UYVY
        {sizeof(DDPIXELFORMAT), DDPF_FOURCC, mmioFOURCC('Y','U','Y','2'), 16, 0, 0, 0, 0}, // YUY2 is unsupported color-space here
        {0}};

    // Setup the overlay surface's attributes in the surface descriptor
    memset(&ddsd, 0, sizeof(ddsd));
    ddsd.dwSize = sizeof(ddsd);
    ddsd.ddsCaps.dwCaps = DDSCAPS_OVERLAY | g_DDCaps.ddsCaps.dwCaps&DDSCAPS_VIDEOMEMORY;
    ddsd.dwFlags = DDSD_CAPS | DDSD_HEIGHT | DDSD_WIDTH | DDSD_PIXELFORMAT;
    ddsd.dwBackBufferCount = 0;
    ddsd.dwWidth = g_sizex;
    ddsd.dwHeight = g_sizey;
    for(int format = 0; ddpfOverlayFormats[format].dwSize; format++) {
        ddsd.ddpfPixelFormat = ddpfOverlayFormats[format];
        // Attempt to create the surface with theses settings
        hRet = g_pDD->CreateSurface(&ddsd, &g_pDDSOverlay, NULL);
        if(hRet == DD_OK) break;
    }
    if (hRet != DD_OK)
        return DisplayError("Can't create appropriate overlay surface", hRet);
    return true;
}

inline void mouse(int k, LPARAM lParam)
{
    int x = (int)LOWORD(lParam), y = (int)HIWORD(lParam);
    g_video->on_mouse( x*g_sizex/(g_rcDst.right - g_rcDst.left),
        y*g_sizey/(g_rcDst.bottom - g_rcDst.top), k);
}

LRESULT CALLBACK InternalWndProc(HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
    PAINTSTRUCT                 ps;         // Structure for the paint message
    POINT                       p = {0, 0}; // Translation point for the window's client region
    HRESULT                     hRet;

    switch (iMsg)
    {
        case WM_MOVE:
            // Make sure we're not moving to be minimized - because otherwise
            // our ratio varialbes (g_dwXRatio and g_dwYRatio) will end up
            // being 0, and once we hit CheckBoundries it divides by 0.
            if (!IsIconic(hwnd))
            {
                g_rcSrc.left = 0;
                g_rcSrc.right = g_sizex;
                g_rcSrc.top = 0;
                g_rcSrc.bottom = g_sizey;
                GetClientRect(hwnd, &g_rcDst);
                g_dwXRatio = (g_rcDst.right - g_rcDst.left) * 1000 /
                             (g_rcSrc.right - g_rcSrc.left);
                g_dwYRatio = (g_rcDst.bottom - g_rcDst.top) * 1000 /
                             (g_rcSrc.bottom - g_rcSrc.top);
                ClientToScreen(hwnd, &p);
                g_rcDst.left = p.x;
                g_rcDst.top = p.y;
                g_rcDst.bottom += p.y;
                g_rcDst.right += p.x;
                CheckBoundries();
            }
            else
                // Else, hide the overlay... just in case we can't do
                // destination color keying, this will pull the overlay
                // off of the screen for the user.
                if (g_pDDSOverlay && g_pDDSPrimary)
                    g_pDDSOverlay->UpdateOverlay(NULL, g_pDDSPrimary, NULL, DDOVER_HIDE, NULL);
            // Check to make sure our window exists before we tell it to
            // repaint. This will fail the first time (while the window is being created).
            if (hwnd)
            {
                InvalidateRect(hwnd, NULL, FALSE);
                UpdateWindow(hwnd);
            }
            return 0L;

        case WM_SIZE:
            // Another check for the minimization action.  This check is
            // quicker though...
            if (wParam != SIZE_MINIMIZED)
            {
                GetClientRect(hwnd, &g_rcDst);
                ClientToScreen(hwnd, &p);
                g_rcDst.left = p.x;
                g_rcDst.top = p.y;
                g_rcDst.bottom += p.y;
                g_rcDst.right += p.x;
                g_rcSrc.left = 0;
                g_rcSrc.right = g_sizex;
                g_rcSrc.top = 0;
                g_rcSrc.bottom = g_sizey;
                // Here we multiply by 1000 to preserve 3 decimal places in the
                // division opperation (we picked 1000 to be on the same order
                // of magnitude as the stretch factor for easier comparisons)
                g_dwXRatio = (g_rcDst.right - g_rcDst.left) * 1000 /
                             (g_rcSrc.right - g_rcSrc.left);
                g_dwYRatio = (g_rcDst.bottom - g_rcDst.top) * 1000 /
                             (g_rcSrc.bottom - g_rcSrc.top);
                CheckBoundries();
            }
            return 0L;

        case WM_PAINT:
            BeginPaint(hwnd, &ps);
            // Check the primary surface to see if it's lost - if so you can
            // pretty much bet that the other surfaces are also lost - thus
            // restore EVERYTHING!  If we got our surfaces stolen by a full
            // screen app - then we'll destroy our primary - and won't be able
            // to initialize it again. When we get our next paint message (the
            // full screen app closed for example) we'll want to try to reinit
            // the surfaces again - that's why there is a check for
            // g_pDDSPrimary == NULL.  The other option, is that our program
            // went through this process, could init the primary again, but it
            // couldn't init the overlay, that's why there's a third check for
            // g_pDDSOverlay == NULL.  Make sure that the check for
            // !g_pDDSPrimary is BEFORE the IsLost call - that way if the
            // pointer is NULL (ie. !g_pDDSPrimary is TRUE) - the compiler
            // won't try to evaluate the IsLost function (which, since the
            // g_pDDSPrimary surface is NULL, would be bad...).
            if (!g_pDDSPrimary || (g_pDDSPrimary->IsLost() != DD_OK) ||
                (g_pDDSOverlay == NULL))
            {
                DestroyOverlay();
                DestroyPrimary();
                if (DDPrimaryInit())
                    if (DDOverlayInit())
                        if (!DrawOverlay())
                            DestroyOverlay();
            }
            // UpdateOverlay is how we put the overlay on the screen.
            if (g_pDDSOverlay && g_pDDSPrimary && g_video->updating)
            {
                hRet = g_pDDSOverlay->UpdateOverlay(&g_rcSrc, g_pDDSPrimary,
                                                    &g_rcDst, g_OverlayFlags,
                                                    &g_OverlayFX);
#ifdef _DEBUG
                if(hRet != DD_OK) DisplayError("Can't update overlay", hRet);
#endif
            }
            EndPaint(hwnd, &ps);
            return 0L;

        // process mouse and keyboard events
        case WM_LBUTTONDOWN:    mouse(1, lParam); break;
        case WM_LBUTTONUP:      mouse(-1, lParam); break;
        case WM_RBUTTONDOWN:    mouse(2, lParam); break;
        case WM_RBUTTONUP:      mouse(-2, lParam); break;
        case WM_MBUTTONDOWN:    mouse(3, lParam); break;
        case WM_MBUTTONUP:      mouse(-3, lParam); break;
        case WM_CHAR:           g_video->on_key(wParam); break;

        case WM_DISPLAYCHANGE:  return 0L;

        case WM_DESTROY:
            // Now, shut down the window...
            PostQuitMessage(0);
            return 0L;
    }
    return g_pUserProc? g_pUserProc(hwnd, iMsg, wParam, lParam) : DefWindowProc(hwnd, iMsg, wParam, lParam);
}

DWORD WINAPI thread_vsync(LPVOID lpParameter)
{
    BOOL vblank = false;
    while(g_video && g_video->running) {
        while(!vblank && g_video && g_video->running) {
            YIELD_TO_THREAD();
            LPDIRECTDRAW7 pDD(g_pDD);
            if(pDD) pDD->GetVerticalBlankStatus(&vblank);
        }
        LPDIRECTDRAWSURFACE7 pDDSOverlay(g_pDDSOverlay);
        if(pDDSOverlay) pDDSOverlay->UpdateOverlay(&g_rcSrc, g_pDDSPrimary, &g_rcDst, g_OverlayFlags | DDOVER_REFRESHALL, &g_OverlayFX);
        do {
            Sleep(1);
            LPDIRECTDRAW7 pDD(g_pDD);
            if(pDD) pDD->GetVerticalBlankStatus(&vblank);
        } while(vblank && g_video && g_video->running);
        while(g_video && !g_video->updating && g_video->running) Sleep(10);
    }
    return 0;
}

///////////////////////////////////////////// public methods of video class ///////////////////////

inline void mask2bits(unsigned int mask, color_t &save, depth_t &shift)
{
    save  = mask; if(!mask) { shift = 8; return; }
    shift = 0; while(!(mask&1)) ++shift, mask >>= 1;
    int bits = 0; while(mask&1) ++bits,  mask >>= 1;
    shift += bits - 8;
}

bool video::init_window(int sizex, int sizey)
{
    assert(win_hInstance != 0);
    g_sizex = sizex; g_sizey = sizey;
    if( !WinInit(win_hInstance, win_iCmdShow, gWndClass, title, false) )
        return DisplayError("Unable to initialize the program's window.");
    running = true;
    if( !DDInit() ) {
        DestroyDDraw();
        goto fail;
    }
    if( !DDOverlayInit() || !DrawOverlay() ) {
        DestroyOverlay();
        DestroyDDraw();
        goto fail;
    }
    DDPIXELFORMAT PixelFormat; memset(&PixelFormat, 0, sizeof(PixelFormat)); PixelFormat.dwSize = sizeof(PixelFormat);
    g_pDDSOverlay->GetPixelFormat(&PixelFormat);
    mask2bits(PixelFormat.dwRBitMask, red_mask, red_shift);
    mask2bits(PixelFormat.dwGBitMask, green_mask, green_shift);
    mask2bits(PixelFormat.dwBBitMask, blue_mask, blue_shift);
    if(PixelFormat.dwFlags == DDPF_RGB)
         depth = depth_t(PixelFormat.dwRGBBitCount);
    else depth = -depth_t(PixelFormat.dwFourCC);
    for(int i = 0, e = sizex * sizey * PixelFormat.dwRGBBitCount / 32, c = get_color(0, 0, 0); i < e; i++)
        g_pImg[i] = c; // clear surface
    ShowWindow(g_hAppWnd, SW_SHOW);
    g_hVSync = CreateThread (
        NULL,          // LPSECURITY_ATTRIBUTES security_attrs
        0,             // SIZE_T stacksize
        (LPTHREAD_START_ROUTINE) thread_vsync,
        this,               // argument
        0, 0);
    SetPriorityClass(g_hVSync, IDLE_PRIORITY_CLASS); // questionable
    return true;
fail:
    g_pImg = new unsigned int[g_sizex * g_sizey];
    return false;
}

void video::terminate()
{
    running = false;
    DestroyOverlay();
    if(WaitForSingleObject(g_hVSync, 100) == WAIT_TIMEOUT) TerminateThread(g_hVSync, 0);
    CloseHandle(g_hVSync);
    DestroyDDraw();
    if(g_pImg) delete[] g_pImg;
    g_pImg = 0; g_video = 0;
}
//////////// drawing area constructor & destructor /////////////

drawing_area::drawing_area(int x, int y, int sizex, int sizey)
: start_x(x), start_y(y), size_x(sizex), size_y(sizey), pixel_depth(g_video->depth),
    base_index(y*g_sizex + x), max_index(g_sizex*g_sizey), index_stride(g_sizex), ptr32(g_pImg)
{
    assert(ptr32); assert(x < g_sizex); assert(y < g_sizey);
    assert(x+sizex <= g_sizex); assert(y+sizey <= g_sizey);

    index = base_index; // current index
}

void drawing_area::update()
{
}

#endif //_DXSDK_PRODUCT_MAJOR >= 9
