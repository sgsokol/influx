#!/usr/bin/env python
"""influx_idyn.py produces dynamic visualisation of isotopomer labeling
propagation in a metabolic network

usage: influx_idyn.py [-h|--help] [--DEBUG] [network[.ftbl]]
"""
#Todo:
#+transparence
#+zooming
#+OnSize
#+video controls
#-arrows
#-data table (plot?)

# 2009-03-24 sokol
# Copyright 2009, MetaSys/INRA, France

# This file content:
# -imports
# -custom classess
# -config constants
# -global vars
# -call backs defs
# -working functions
# -line arguments parse
# -GUI layout (from *_lay.kvh)

## imports
import sys;
import os;
import time;
import getopt;
import re;
import numpy as np;
import array;
import math;
from threading import Thread;

import wx;
import wx.html;
import wx.grid as wxg;
import wx.lib.ogl as wxo;

import webbrowser;
import xml.dom.minidom as xd;

#sys.path.append('/home/sokol/dev/python/pytools');
#sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
dirx=os.path.dirname(sys.argv[0]) or os.getcwd();
sys.path.append(dirx);
import tools_ssg as tls;
#from C13_ftbl import *;
import C13_ftbl as c13;

## custom classes
class data2tab(wxg.PyGridTableBase):
    def __init__(self, data, cnames=[]):
        wxg.PyGridTableBase.__init__(self);
        self.data=data;
        self.cnames=cnames;
    def GetNumberRows(self):
        return(self.data.shape[0]);
    def GetNumberCols(self):
        return(self.data.shape[1]);
    def IsEmptyCell(self, row, col):
        return False
    def GetValue(self, row, col):
        return("%.3f"%self.data[row, col]);
        #return(self.data[row, col]);
    def SetValue(self, row, col, value):
        pass;
    def GetColLabelValue(self, col):
        return self.cnames[col] if self.cnames else None;

class fitGrid(Thread):
    def __init__(self, grid):
        Thread.__init__(self);
        self.grid=grid;
    def run(self):
        self.grid.Fit();
        parent=self.grid.GetParent();
        parent.SetVirtualSize(self.grid.GetSize());
        #parent.Refresh();
        #parent.Enable(True);

## config constants
# program name
me=os.path.basename(sys.argv[0]);
# border width in pixels for widgets
wd_border=3;
tsize=(24,24); # tollbar size
ID_PLAY=10; # id for play button
ID_PAUSE=20; # id for pause button
# message in welcome tab
welc_text="""
<html>
<h1>Welcome to %(me)s!</h1>

<p>
%(me)s is part of influx_sim package designed for dynamic presentation of
results in isotopomer non stationnary simulations.
</p>

<h3>Usage.</h3>
<ol>
<li>Open an FTBL file: File > Open > ...(choose some file, let call it network.ftbl)

(FTBL is a format defined in <a href="http://www.uni-siegen.de/fb11/simtec/software/13cflux/">
13cflux</a> software written by W. Wiechert.)



There are some name conventions to respect for file names. Thus,
for a file network.ftbl it must be in the same directory the following files:
<ul>
<li>network.xgmml (network graph produced by <a href="http://www.cytoscape.org">cytoscape</a>)</li>
<li>network_ipl.kvh (data to plot produced by <tt>influx_i.sh network.ftbl</tt>) </li>
<li>network_idyn.kvh (your custom settings for the network saved
from previous session)</li>
</ul>
</li>
<li>Data to plot are presented in <i>Data</i> tab.
Check if it is the data you mean to plot.
</li>
<li>Network layout is presented in <i>Network</i> tab according to xgmml
content. Use video controls to launch the dynamic presentation
which will use colour codes on metabolite nodes to reflect the labelling
propagation in the network
</li>
<li>Use <i>Settings</i> tab to customize colour codes, video speed and so on.
</li>
</ol>
<p>
For legal information see Help > About
</p>
</html>
""" % {"me": me};
welc_text=re.sub("\n\n", "<br>\n", welc_text);

## global vars
network="";
pnetw="";
cnames=[];
data=np.array([]);
settings={
    "mrg_x": 10,
    "mrg_y": 10,
    "radius": -0.2, # negative radius to make it proportional to the lowest size of rectangle
    "fontPointSize": 10,
    "tw2tm": 1., # coefficient for time conversion from wallclock to modeled time
    "colour_bar_bold": "#346CDF",
    "frame_period": 100, # milliseconds between two updates of the picture
    "zoom": 1, # zooming factor
};
#pnodes={}; # node plot attributes
pmetabs=set(); # set of plotted metabolites
dom=None;
#xnodes=None;
netw_img=None; # static image of the network
graph={}; # dictionary describing metab graph (node, edge and so on)
ti_wcs=None; # wall clocktime at the model's start
ti_ms=None; # model's start time (after a pause it may be non zero)
base_rect=None; # list of (x,y,w) coord of left lower corner for each bar
dind=None; # dict of column indexes for each metabolite

## call back functions
def OnExit(evt):
    """
    This is executed when the user clicks the 'Exit' option
    under the 'File' menu or close the window.  We ask the user if he *really*
    want to exit, then close everything down if he does.
    """
    #for a in dir(evt):
    #    if a[:3]=="Get":
    #        print("evt."+a+"="+eval("str(evt."+a+"())"));
    win=evt.GetEventObject();
    if "TopLevelParent" in dir(win):
        win=win.TopLevelParent;
    elif "InvokingWindow" in dir(win):
        win=win.InvokingWindow;
    else:
        dlg=wx.MessageDialog(None, "Oups. Should not be there.\nUknown event type in OnExit()", 'Error', wx.OK | 
            wx.ICON_ERROR);
        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy();
            sys.exit(1);
    dlg = wx.MessageDialog(win, 'Exit such a beautifull program?', 'I Need To Know!',
                          wx.YES_NO | wx.ICON_QUESTION);
    if dlg.ShowModal() == wx.ID_YES:
        dlg.Destroy();
        win.Destroy();
    else:
        dlg.Destroy();
def OnOpen(evt):
    """
    This is executed when the user clicks the 'Open' option
    under the 'File' menu.  We ask the user to choose an ftbl file.
    """
    win=evt.GetEventObject();
    win=win.InvokingWindow;
    dlg = wx.FileDialog(win, defaultDir=os.getcwd(), wildcard="FTBL files (*.ftbl)|*.ftbl",
        style=wx.OPEN);
    if dlg.ShowModal() == wx.ID_OK:
        #print "selected file="+dlg.GetPath();
        # proceed the ftbl
        get_proj(dlg.GetPath(), win);
    else:
        dlg.Destroy();
def OnSave(evt):
    """
    This is executed when the user clicks the 'Save settings' option
    under the 'File' menu. Settings are stored in network_idyn.kvh.
    """
    if not network:
        # network is not yet choosed
        dlg=wx.MessageDialog(None, "No network yet chosen.\nChoose a FTBL file first.", 'Error', wx.OK | 
            wx.ICON_ERROR);
        dlg.ShowModal();
        dlg.Destroy();
        return;
    # file name
    fn=network+"_idyn.kvh";
    fpath=os.path.join(pnetw,fn);
    if os.path.exists(fpath):
        # ask permission to overwrite
        dlg=wx.MessageDialog(None, "Overwrite "+fpath+"?", "Confirm", 
                wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION);
        if dlg.ShowModal() == "Yes":
            # write sttings
            dict2kvh(settings, fpath);
        else:
            dlg.Destroy();
    else:
        # write sttings
        dict2kvh(settings, fpath);

def OnPaint(evt):
    #print("pe=", evt);
    #sw_netw.SetBackgroundStyle(wx.BG_STYLE_CUSTOM);
    if netw_img:
        pbuf=netw_img.ConvertToBitmap();
        dc = wx.BufferedPaintDC(pan_netw, pbuf);
        #print("pos=", evt.GetEventObject().GetViewStart());
        #vs=evt.GetEventObject().GetViewStart()
        #dc.DrawBitmap(pbuf, -vs[0]*20, -vs[1]*20);

def OnLinkClicked(evt):
    #print(str(dir(evt)));
    #for o in dir(evt.GetLinkInfo()):
    #    print(o+"="+str(eval("evt.GetLinkInfo()."+o)));
    webbrowser.open_new_tab(evt.GetLinkInfo().Href);

def OnTimer(evt):
    """plot the current state of labeled matabs
    if evt==None, set timers to zero.
    """
    global ti_wcs, ti_ms, netw_img;
    if ti_wcs is None:
        ti_wcs=time.time();
    if ti_ms is None:
        ti_ms=data[0,0]; # starting value for modeled time
    ti_wcc=time.time();
    ti_m=ti_ms+(ti_wcc-ti_wcs)*settings["tw2tm"];
    ftime=s2ftime(ti_m);
    pbuf=netw_img.ConvertToBitmap();
    dc=wx.BufferedPaintDC(pan_netw, pbuf);
    #dc.DrawLabel(ftime, (ti_wcc-ti_wcs, ti_wcc-ti_wcs, 100, 100));
    #print("ftime=", ftime);
    labti_m.SetLabel(ftime);
    drawBars(dc, ti_m);
    if ti_m > data[-1,0]:
        # we are at the end of modeled time
        timer.Stop();
        ti_ms=data[0,0];
        netw_img=pbuf.ConvertToImage();
        pan_netw.Refresh();
        play.Enable(True);
        pause.Enable(False);

def OnPlay(evt):
    """launch the animation
    """
    global ti_wcs, ti_ms;
    ti_wcs=time.time();
    if ti_ms is None:
        ti_ms=data[0,0]; # starting value for modeled time
    timer.Start(settings["frame_period"]);
    play.Enable(False);
    pause.Enable(True);

def OnPause(evt):
    "pause the animation"
    global ti_ms;
    timer.Stop();
    ti_wcc=time.time();
    ti_ms=ti_ms+(ti_wcc-ti_wcs)*settings["tw2tm"];
    play.Enable(True);
    pause.Enable(False);

def OnSize(evt):
    "main window is resized"
    win=evt.GetEventObject();
    sz=evt.GetSize();
    #print("ons=", win, sz);
    win.SetSize(sz);
    evt.Skip();
    return;

def OnSlider(evt):
    global settings, ti_ms, ti_wcs;
    if evt.GetEventObject().Id == speed.Id:
        v=speed.GetValue();
        f=4**((v-50.)/50.);
        settings["tw2tm"]=f;
        lf="%.1f"%f;
        labspeed.SetLabel("" if lf=="1.0" else "x"+lf if f > 1. else "/%.1f"%(1./f));
    if evt.GetEventObject().Id == vpos.Id:
        v=vpos.GetValue();
        #print("v=",v);
        # set the modeled time at this position
        ti_wcs=time.time();
        ti_ms=v;
        labti_m.SetLabel(s2ftime(ti_ms));
        if play.Enabled:
            #print("try to draw");
            pbuf=netw_img.ConvertToBitmap();
            drawBars(wx.BufferedPaintDC(pan_netw, pbuf), ti_ms);

def ToDo(evt):
    """
    A general purpose "we'll do it later" dialog box
    """
    win=evt.GetEventObject().GetTopLevelParent();
    dlg = wx.MessageDialog(win, "Not Yet Implimented! evt="+str(dir(evt)), "ToDo",
                         wx.OK | wx.ICON_INFORMATION);
    dlg.ShowModal();
    dlg.Destroy();

## working functions
def get_proj(fn, mf):
    """get_proj(fn)
    parse files associated to ftbl given in fn
    """
    global network, pnetw, cnames, data, dom, netan, pmetabs,\
        netw_img, dind, graph, ti_ms;
    if not fn or fn[-5:] != ".ftbl":
        return;
    # parse ftbl to netan
    mainframe.SetStatusText("Parsing "+fn)
    #print("ftbl=",c13.ftbl_parse(fn));
    netan=c13.ftbl_netan(c13.ftbl_parse(fn));
    
    # set file names
    network=os.path.basename(fn)[:-5];
    pnetw=os.path.dirname(fn);
    
    fdata=os.path.join(pnetw, network+"_ipl.kvh");
    fx=os.path.join(pnetw, network+".xgmml");
    #print ("fdata=", fdata);
    # get data to plot
    (cnames,data)=ipl2data(fdata);
    cnames[0]="time";
    ti_ms=data[0,0];
    # set min max video position
    vpos.SetMin(int(ti_ms));
    vpos.SetMax(int(data[-1,0]));
    # fill the set of metabolite names
    pmetabs=set(re.split("[:#+]", s)[0] for s in cnames[1:]);
    pmetabs.update(netan["input"]);
    # store indeces of data columns for each metabolite
    dind={};
    for (ic,nm) in enumerate(cnames):
        m=re.split("[:#+]", nm)[0];
        dind[m]=dind.get(m, []);
        dind[m].append(ic);
    #print("dind=", dind);
    
    # parse xml
    dom=xd.parse(fx);
    # parse graph
    graph=dom2gr(dom);
    #print("dom=", dom);
    
    ## create graph layout
    # create static buffer for the network graph
    mrg_x=settings["mrg_x"];
    mrg_y=settings["mrg_y"];
    (w,h)=(graph["box"][-2:]);
    (w,h)=(int(w+2*mrg_x), int(h+2*mrg_y));
    #pbuf=wx.EmptyBitmap(w+2*mrg_x, h+2*mrg_y);
    bpp = 4  # bytes per pixel
    bytes = array.array("B", [128] * (w*h*bpp));
    pbuf=wx.BitmapFromBufferRGBA(w, h, bytes);
    # connect the buffer to DC
    dc=wx.BufferedPaintDC(pan_netw, pbuf);
    # draw to buffer via DC
    drawGraph(graph, dc);
    # draw the starting labeling
    drawBars(dc, ti_ms);
    # store the static image
    netw_img=pbuf.ConvertToImage();
    play.Enable();
    i,p=lab2ip(nb, "Network");
    nb.SetSelection(i); # show directly network tab
    #print(type(pbuf));
    
    ## create and fill the data table
    if "grid" in dir(sw_data):
        # destroy the previous grid
        sw_data.grid.Destroy();
    # create new one
    grid=wxg.Grid(sw_data, wx.ID_ANY);
    sw_data.grid=grid;
    grid.table=data2tab(data, cnames);
    #print("gt=", g.table);
    grid.SetTable(grid.table, True);
    # before going to thread, disable access to data tab
    # if not, it hangs
    #i,p=lab2ip(nb, "Data");
    #p.Enable(False); # the panel is enabled at the end of thread
    th=fitGrid(grid);
    th.start();

def ipl2data(fn):
    fp=open(fn, "r");
    # get column names
    cnames=fp.readline().strip().split("\t");
    #print cnames;
    #print fp.tell();
    data=np.loadtxt(fp);
    #print data;
    return (cnames, data);

def drawGraph(graph, dc):
    mrg_x=settings["mrg_x"];
    mrg_y=settings["mrg_y"];
    dc.SetPen(wx.Pen("BLACK",1));
    (w,h)=(graph["box"][-2:]);
    pan_netw.SetSize(((w+2*mrg_x)*settings["zoom"], (h+2*mrg_y)*settings["zoom"]));
    sw_netw.SetVirtualSize(pan_netw.GetSize());
    dc.Clear();
    #dc.SetLogicalOrigin(-mrg_x,-mrg_y);
    # set font size
    font=dc.GetFont();
    font.SetPointSize(int(settings["fontPointSize"]));
    dc.SetFont(font);
    dc.SetUserScale(settings["zoom"], settings["zoom"]);
    #canvas=wxo.ShapeCanvas(sw_netw);
    # draw edges
    n=graph["node"];
    idl=graph["id2lab"];
    z=settings["zoom"];
    for e in graph["edge"]:
        s=n[idl[e["source"]]];
        t=n[idl[e["target"]]];
        b=e.get("edgeBend");
        lp=[(s["x"]+s["w"]/2, s["y"]+s["h"]/2),
                (t["x"]+t["w"]/2, t["y"]+t["h"]/2)]
        if b:
            # draw bended line
            lp=lp[0:1]+\
                [(i["x"], i["y"]) for i in b]+\
                lp[-1:];
        #lp=[(x*z, y*z) for (x,y) in lp];
        #InsertLineControlPoint(self, dc=None, point=None)
        dc.DrawSpline(lp);
    # draw strait rectangles
    #print("nodes (x,y,w,h)=", [(pn["x"], pn["y"], pn["w"], pn["h"]) for pn in pnodes if pn["shape"]=="RECTANGLE"]);
    dc.DrawRectangleList([(pn["x"], pn["y"], pn["w"], pn["h"])
        for pn in graph["node"].values() if pn["label"] in pmetabs or pn["gr"]["type"]=="RECTANGLE"]);
    # draw nodes (rounded rectangles, rombs, etc.)
    # and their labels
    pen=dc.GetPen();
    for pn in graph["node"].values():
        if pn["label"] not in pmetabs:
            if pn["gr"]["type"]=="ROUNDED_RECTANGLE":
                dc.SetBrush(wx.Brush(pn["gr"]["fill"]));
                dc.DrawRoundedRectangle(pn["x"], pn["y"], pn["w"], pn["h"], settings["radius"]);
            elif pn["gr"]["type"]=="DIAMOND":
                dc.SetBrush(wx.Brush(pn["gr"]["fill"]));
                #dc.DrawRectangle(pn["x"], pn["y"], pn["w"], pn["h"]);
                dc.DrawPolygon(dia2points(pn["x"], pn["y"], pn["w"], pn["h"]));
        else:
            # inputs are in green
            dc.SetBrush(wx.Brush(wx.WHITE));
            if pn["label"] in netan["input"]:
                dc.SetPen(wx.Pen(wx.GREEN));
                dc.DrawRectangle(pn["x"], pn["y"], pn["w"], pn["h"]);
                #r=wxo.DrawnShape.DrawRectangle((pn["x"], pn["y"], pn["w"], pn["h"]));
                dc.SetPen(pen);
        dc.DrawLabel(pn["label"],
            (pn["x"], pn["y"], pn["w"], pn["h"]), wx.ALIGN_CENTER);

def dom2gr(dom):
    """dom2gr(dom)-> dict
    extract nodes and edges from dom structure (xgmml)
    """
    mrg_x=settings["mrg_x"];
    mrg_y=settings["mrg_y"];
    res={};
    xnodes=dom.getElementsByTagName("node");
    res["node"]=[el2dict(n) for n in xnodes];
    # id to label translation
    res["id2lab"]=dict((val["id"], val["label"]) for val in res["node"]);
    # node list to dict with indexing by node label
    res["node"]=dict((val["label"],val) for val in res["node"]);
    
    #print("xn=", dir(xnodes));
    # get graph node limits
    minx=sys.maxint;
    miny=sys.maxint;
    w=0;
    h=0;
    for n in res["node"].values():
        #print("n=", n);
        n["x"]=float(n["gr"]["x"]);
        n["y"]=float(n["gr"]["y"]);
        n["w"]=float(n["gr"]["w"]);
        n["h"]=float(n["gr"]["h"]);
        minx=min(minx, n["x"]);
        miny=min(miny, n["y"]);
        w=max(w,n["x"]+n["w"]);
        h=max(h,n["y"]+n["h"]);
    
    (w,h)=(w-minx, h-miny);
    dx=-minx+mrg_x;
    dy=-miny+mrg_y;
    res["box"]=(0, 0, w, h);
    #print("gr w,h=", w, h);
    # get graph edges
    res["edge"]=[el2dict(e) for e in dom.getElementsByTagName("edge")];
    # translate all nodes and edges
    for n in res["node"].values():
        n["x"]+=dx;
        n["y"]+=dy;
    for e in res["edge"]:
        b=e.get("edgeBend", []);
        for item in b:
            item["x"]=float(item["x"])+dx;
            item["y"]=float(item["y"])+dy;
    return(res);

def el2dict(e):
    res={};
    if e.hasAttributes():
        res.update(e.attributes.items());
    gr=e.getElementsByTagName("graphics");
    if gr:
        res["gr"]=el2dict(gr[0]);
    ats=e.getElementsByTagName("att");
    if ats:
        res["att"]={};
    for a in ats:
        d=att2dict(a);
        if a.hasChildNodes():
            res[d.get("name")]=[att2dict(suba) for suba in a.getElementsByTagName("att")];
        else:
            res["att"][d.get("name")]=d.get("value");
    return(res);

def att2dict(a):
    """att2dict(a) -> dict
    transofrm <a x=ax y=ay/> in dict("x": "ax", "y": "ay");
    """
    return dict(a.attributes.items());

def dia2points(x, y, w, h):
    """dia2points(x, y, w, h)-> tuple of 4 couples (px,py)
    diamond is inscribed in the rectangle (x,y,w,h)
    x,y is the top left point of the rectangle
    w,h are the width and heght of the rectangle
    """
    #z=settings["zoom"];
    w*=0.5;
    h*=0.5;
    #x*=z;
    #y*=z;
    return ((x, y+h), (x+w,y), (x+w+w, y+h), (x+w, y+h+h));

def s2ftime(s=0.):
    """s2ftime(s=0) -> String
    Format second number as hh:mm:ss.cc
    """
    si=int(s);
    cc=round(100*(s-si), 0);
    s=si;
    ss=s%60;
    s//=60;
    mm=s%60;
    s//=60;
    hh=s;
    return("%02d:%02d:%02d.%02d"%(hh,mm,ss,cc));

def drawBars(dc, ti_m):
    """drawBars(dc, ti_m) -> None
    Draw rectangle bars in the metabolite nodes representing
    labeling propagation.
    """
    global base_rect;
    mrg_x=settings["mrg_x"];
    mrg_y=settings["mrg_y"];
    #dc.SetLogicalOrigin(-mrg_x,-mrg_y);
    dc.SetUserScale(settings["zoom"], settings["zoom"]);
    nds=graph["node"];
    if not base_rect:
        # prepare list of tuples (x,y,w) of lower left corner of each rectangle
        base_rect=[];
        for m in pmetabs:
            if m not in dind:
                continue; # no data for this metabolite, it must be an input
            nd=nds[m];
            py=nd["y"]+nd["h"];
            n=len(dind[m]);
            w=nd["w"]/n;
            pxs=[round(nd["x"]+i*w) for i in xrange(n)];
            #print("m=", m, "pxs=", pxs);
            for i in xrange(n-1):
                base_rect.append((pxs[i],py,pxs[i+1]-pxs[i], m));
            # treat the last bar
            base_rect.append((pxs[n-1], py, pxs[0]+nd["w"]-pxs[n-1], m));
    #print("br=", base_rect);
    #print("n=", [(nds[l]["x"], nds[l]["y"], nds[l]["w"], l) for l in pmetabs]);
    #sys.exit(1);
    # find the max row index in data s.t ti(ir)<=tm
    # and fraction f of ir+1 s.t. 0<=f<=1, f=0 => 100% ir, f=1 => 100% of ir+1
    (nr,nc)=data.shape;
    if ti_m >= data[-1,0]:
        ir=nr-2;
        f=1.;
    elif ti_m <= data[0,0]:
        ir=0.;
        f=0.;
    else:
        for ir in xrange(nr):
            if data[ir,0]>ti_m:
                ir=ir-1;
                f=(ti_m-data[ir,0])/(data[ir+1,0]-data[ir,0]);
                break;
    # prepare the heights of all bars
    # the bars ar proportional to the faction of labeled component
    #print("ti_m=", ti_m, "ir=", ir, "f=", f);
    hs=[];
    for m in pmetabs:
        if m not in dind:
            continue; # no data for this metabolite, it must be an input
        nd=nds[m];
        n=len(dind[m]);
        h=nd["h"];
        for i in xrange(n):
            ic=dind[m][i];
            v=data[ir,ic]*(1.-f)+data[ir+1,ic]*f;
            hs.append(round(h*v));
    # gather all ractangles and draw them
    rect=[(b[0],b[1]-hs[i],b[2],hs[i]) for (i,b) in enumerate(base_rect)];
    dc.SetPen(wx.Pen(settings["colour_bar_bold"], 0));
    brush=wx.Brush(settings["colour_bar_bold"]);
    (r,g,b)=brush.GetColour().Get()
    dc.SetBrush(wx.Brush(wx.Colour(r,g,b,alpha=128)));
    #dc.SetBrush(wx.Brush(settings["colour_bar_bold"], wx.TRANSPARENT|wx.SOLID));
    DrawRectangleListAsBitmap(dc, rect);

def DrawRectangleListAsBitmap(dc, rect):
    """DrawRectangleListAsBitmap(dc, rect)->None
    Create a rectangle bitmap for each rectangle to
    be able to use the transparency.
    """
    # get max w,h
    (maxw,maxh)=(0,0);
    for r in rect:
        (maxw,maxh)=(max(maxw,r[2]), max(maxh,r[3]));
    (maxw,maxh)=(int(maxw),int(maxh));
    # allocate buffer with brush colour values
    colour=dc.GetBrush().GetColour();
    (r,g,b)=colour;
    a=colour.Alpha();
    #print("a=", a);
    buf=array.array("B", [r,g,b,a]*(maxw*maxh));
    for (x,y,w,h) in rect:
        (x,y,w,h)=(int(x),int(y),int(w),int(h));
        #print("x,y,w,h=",x,y,w,h);
        if w*h <= 0:
            continue;
        bmp=wx.BitmapFromBufferRGBA(w, h, buf[:4*w*h]);
        dc.DrawBitmap(bmp, x, y, True);

def lab2ip(nb, lab):
    """lab2i(nb, nm) -> (i,Page) or (None,None)
    get page of a notebook nb by its label lab
    """
    nb_pages=nb.GetPageCount();
    for i in xrange(nb_pages):
        if lab == nb.GetPageText(i):
            return((i,nb.GetPage(i)));
    return (None,None);
## take arguments
#<--skip in interactive session
# get arguments
def usage():
    print(__doc__);
try:
    opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "DEBUG"]);
except getopt.GetoptError, err:
    print str(err);
    usage();
    sys.exit(1);
DEBUG=False;
for o,a in opts:
    if o in ("-h", "--help"):
        usage();
        sys.exit(0);
    elif o=="--DEBUG":
        DEBUG=True;
    else:
        assert False, "unhandled option";
fftbl=args[0] if len(args) else "";
if fftbl and fftbl[-5:] != ".ftbl":
    fftbl+=".ftbl";
if fftbl and not os.path.exists(fftbl):
    sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n");
    sys.exit(1);
if fftbl:
    wd=os.path.dirname(fftbl);
    if wd != "":
        #print("wd='"+wd+"'");
        os.chdir(wd);
    fftbl=os.path.basename(fftbl);

## create GUI
app=wx.App();
code=tls.wxlay2py(tls.kvh2tlist(sys.argv[0][:-3]+"_lay.kvh"));
flname=sys.argv[0][:-3]+"_lay.py";
fp=open(flname, "w");
fp.write(code);
fp.close();
execfile(flname);
if fftbl:
    get_proj(os.path.join(os.getcwd(), fftbl), mainframe);
app.MainLoop();

