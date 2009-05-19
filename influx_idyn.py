#!/usr/bin/env python
"""influx_idyn.py produces dynamic visualisation of isotopomer labeling
propagation in a metabolic network

usage: influx_idyn.py [-h|--help] [--DEBUG] [network[.ftbl]]
"""
#Todo:
#+transparence
#+OnSize
#+video controls
#+data table
#-arrows
#-(plot?)

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
#import wx.lib.ogl as wxo;
from wx.lib.wordwrap import wordwrap;
import wx.lib.floatcanvas.FloatCanvas as wxfc;

import webbrowser;
import xml.dom.minidom as xd;

#sys.path.append('/home/sokol/dev/python/pytools');
#sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
dirx=os.path.dirname(os.path.abspath(os.sys.argv[0]));
#print("dirx=", dirx);
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

class MovingObjectMixin:
    """
    Methods required for a Moving object
    
    """
    def GetOutlinePoints(self):
        BB = self.BoundingBox
        OutlinePoints = np.array( ( (BB[0,0], BB[0,1]),
                                    (BB[0,0], BB[1,1]),
                                    (BB[1,0], BB[1,1]),
                                    (BB[1,0], BB[0,1]),
                                 )
                               )
        return OutlinePoints

class MovingRectangle(wxfc.Rectangle, MovingObjectMixin):
    """
    ScaledBitmap Object that can be moved
    """
    ## All we need to do is is inherit from:
    ##  ScaledBitmap, MovingObjectMixin
    def GetXYWH(self):
        (x,y),(w,h)=self.XY,self.WH;
        return(x,y,w,h);

# workaround for transparency with gcdc
Draw_orig=wxfc.Rectangle._Draw;
def gcdcDraw(self, dc, *args, **kwargs):
    try:
        dc=wx.GCDC(dc);
        #print "gcdc ok";
    except:
        #print "gcdc ko";
        pass;
    Draw_orig(self, dc, *args, **kwargs);
wxfc.Rectangle._Draw=gcdcDraw;

# work around for negative height
SetUpDraw_orig=wxfc.Rectangle.SetUpDraw;
def SetUpDraw_pos(self, *args, **kwargs):
    (XY, WH) = SetUpDraw_orig(self, *args, **kwargs);
    if WH[1] < 0:
       XY[1] += WH[1]-1;
       WH[1] = np.abs(WH[1]);
    return(XY, WH);
wxfc.Rectangle.SetUpDraw=SetUpDraw_pos;

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

<h3>Usage:</h3>
<ol>
<li>Open an FTBL file: File > Open > ...(choose some file, let call it network.ftbl)

(FTBL is a format defined in <a href="http://www.uni-siegen.de/fb11/simtec/software/13cflux/">
13cflux</a> software written by W. Wiechert.)



There are some name conventions to respect for file names. Thus,
for a file network.ftbl it must be in the same directory the following files:
<ul>
<li>network.xgmml (network graph produced by <a href="http://www.cytoscape.org">cytoscape</a>)</li>
<li>network_ipl.kvh (data to plot produced by <tt>influx_i.sh network.ftbl</tt>) </li>
<li>network_idyn.kvh (your optional custom settings for the network saved
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
fp=open(os.path.join(dirx,"licence_en.txt"), "r");
licenseText=fp.read();
fp.close();

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
};
#pnodes={}; # node plot attributes
pmetabs=set(); # set of plotted metabolites
dom=None;
#xnodes=None;
netw_img=None; # static image of the network
graph={}; # dictionary describing metab graph (node, edge and so on)
ti_wcs=None; # wall clocktime at the model's start
ti_ms=None; # model's start time (after a pause it may be non zero)
bars=None; # list of [x,y,w,h, shape] indexed by column number in data array
dind=None; # dict of column indexes for each metabolite
canvas=None; # plotting canvas for network and animations
sett_inp=None; # dictionary of input widgets for settings

## call back functions
def OnExit(evt):
    """
    This is executed when the user clicks the 'Exit' option
    under the 'File' menu or close the window.  We ask the user if he *really*
    want to exit, then close everything down if he does.
    """
    dlg = wx.MessageDialog(None, 'Exit such a beautifull program?', 'I Need To Know!',
                          wx.YES_NO | wx.ICON_QUESTION);
    if dlg.ShowModal() == wx.ID_YES:
        dlg.Destroy();
        mainframe.Destroy();
    else:
        dlg.Destroy();
def OnOpen(evt):
    """
    This is executed when the user clicks the 'Open' option
    under the 'File' menu.  We ask the user to choose an ftbl file.
    """
    #print("wd=", wd);
    dlg = wx.FileDialog(None, defaultDir=wd, wildcard="FTBL files (*.ftbl)|*.ftbl",
        style=wx.OPEN);
    if dlg.ShowModal() == wx.ID_OK:
        #print "selected file="+dlg.GetPath();
        # proceed the ftbl
        get_proj(dlg.GetPath());
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
        rep=dlg.ShowModal();
        #print("rep=", rep);
        if rep == wx.ID_YES:
            # write sttings
            #print("writing "+fpath);
            tls.dict2kvh(settings, fpath);
            mainframe.SetStatusText("Written "+fpath);
        else:
            dlg.Destroy();
    else:
        # write directly
        tls.dict2kvh(settings, fpath);
        mainframe.SetStatusText("Written "+fpath);

def OnAbout(evt):
    "show about dialog"
    win=evt.GetEventObject();
    win=win.InvokingWindow;
    info = wx.AboutDialogInfo();
    #print("sz=", info.GetSizer());
    info.SetIcon(wx.IconFromLocation(wx.IconLocation(os.path.join(dirx, "icons", "metasys.png"))));
    info.Name = "influx_idyn.py";
    info.Version = "1.0.0";
    info.Copyright = "(C) 2009 INRA";
    info.Description = wordwrap(
        "A 'influx_idyn.py' program is a software that shows"
        " an animation presenting labeling propagation in a metabolic"
        " network which is in metabolicaly stationary state."
        " It comes as visualisation tool for results obtained with"
        " influx_i program." ,
        350, wx.ClientDC(win));
    info.WebSite = ("http://metasys.insa-toulouse.fr", "MetaSys home page");
    info.Developers = [ "Serguei SOKOL" ];

    info.License = wordwrap(licenseText, 500, wx.ClientDC(win));

    # Then we call wx.AboutBox giving it that info object
    wx.AboutBox(info);

def OnLinkClicked(evt):
    #print(str(dir(evt)));
    #for o in dir(evt.GetLinkInfo()):
    #    print(o+"="+str(eval("evt.GetLinkInfo()."+o)));
    webbrowser.open_new_tab(evt.GetLinkInfo().Href);

def OnTimer(evt):
    """plot the current state of labeled matabs
    if evt==None, set timers to zero.
    """
    global ti_wcs, ti_ms;
    if ti_wcs is None:
        ti_wcs=time.time();
    if ti_ms is None:
        ti_ms=data[0,0]; # starting value for modeled time
    ti_wcc=time.time();
    #print("ti_wcc=", ti_wcc);
    ti_m=ti_ms+(ti_wcc-ti_wcs)*settings["tw2tm"];
    ftime=s2ftime(ti_m);
    #print("ti_m, ftime=", ti_m, ftime);
    labti_m.SetLabel(ftime);
    vpos.SetValue(round(ti_m));
    drawBars(ti_m);
    #wx.App.Yield(app);
    if ti_m >= data[-1,0]:
        # we are at the end of modeled time
        timer.Stop();
        ti_ms=data[0,0];
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
    #print("period=", settings["frame_period"]);
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
        show_tfactor(f);
    if evt.GetEventObject().Id == vpos.Id:
        v=vpos.GetValue();
        #print("v=",v);
        # set the modeled time at this position
        ti_wcs=time.time();
        ti_ms=v;
        labti_m.SetLabel(s2ftime(ti_ms));
        if play.Enabled:
            #print("try to draw");
            drawBars(ti_ms);

def OnSettings(evt):
    "store the new value in settings"
    global settings;
    #print(evt.GetString(), dir(inp));
    k=evt.GetEventObject().set_field;
    v=evt.GetString();
    try:
        settings[k]=eval(v);
    except:
        settings[k]=v;
    if k=="frame_period" and timer.IsRunning():
        timer.Stop();
        timer.Start(settings["frame_period"]);
    elif k=="tw2tm":
        try:
            settings[k]=float(settings[k]);
            show_tfactor(settings[k]);
        except:
            pass;
    #print(settings);

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
def get_proj(fn):
    """get_proj(fn)
    parse files associated to ftbl given in fn
    """
    global network, pnetw, cnames, data, dom, netan, pmetabs,\
        dind, graph, ti_ms, settings, wd, bars;
    mf=mainframe;
    if not fn or fn[-5:] != ".ftbl":
        return;
    wd=os.path.dirname(os.path.abspath(fn));
    # parse ftbl to netan
    mainframe.SetStatusText("Parsing "+fn)
    #print("ftbl=",c13.ftbl_parse(fn));
    netan=c13.ftbl_netan(c13.ftbl_parse(fn));
    
    # set file names
    network=os.path.basename(fn)[:-5];
    pnetw=os.path.dirname(fn);
    
    fdata=os.path.join(pnetw, network+"_ipl.kvh");
    fx=os.path.join(pnetw, network+".xgmml");
    fs=os.path.join(pnetw, network+"_idyn.kvh");
    
    # get settings if exist
    if os.path.exists(fs):
        d=tls.kvh2dict(fs);
        for (k,v) in d.iteritems():
            if k not in settings:
                continue;
            try:
                settings[k]=eval(v);
            except:
                settings[k]=v;
    # set slider of speed
    f=settings["tw2tm"];
    show_tfactor(f);
    
    #print("s=", settings);
    
    #print ("fdata=", fdata);
    # get data to plot
    try:
        (cnames,data)=ipl2data(fdata);
    except:
        return;
    # remove precedent project
    bars=None;
    cnames[0]="time";
    ti_ms=data[0,0];
    # set min max video position
    vpos.SetMin(int(ti_ms));
    vpos.SetMax(int(data[-1,0]));
    # fill the set of metabolite names
    pmetabs=set(re.split("[:#+]", s)[0] for s in cnames[1:]);
    pmetabs.update(netan["input"]);
    # store indices of data columns for each metabolite
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
    
    # draw the graph
    drawGraph(graph);
    drawBars(ti_ms);
    
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
    #th=fitGrid(grid);
    #th.start();
    grid.Fit();
    sw_data.SetVirtualSize(grid.GetSize());

def ipl2data(fn):
    try:
        fp=open(fn, "r");
    except:
        dlg=wx.MessageDialog(None, "File "+fn+" could not be read\nChoose another FTBL file or run influx_i before", "Error", wx.OK | wx.ICON_ERROR);
        if dlg.ShowModal() == wx.ID_OK:
            dlg.Destroy();
            return;
    # get column names
    cnames=fp.readline().strip().split("\t");
    #print cnames;
    #print fp.tell();
    data=np.loadtxt(fp);
    #print data;
    return (cnames, data);

def drawGraph(graph):
    global canvas;
    (w,h)=(graph["box"][-2:]);
    mrg_x=settings["mrg_x"];
    mrg_y=settings["mrg_y"];
    w=(w+2*mrg_x);
    h=(h+2*mrg_y);
    # Create new Canvas
    canvas = wxfc.FloatCanvas(pan_netw, wx.ID_ANY,(w,h),
        ProjectionFun = None,
        Debug = 0,
        BackgroundColor = "White",
        );
    #print("w,mrg_x=", (w+2*mrg_x));
    pan_netw.SetSize((w,h));
    sw_netw.SetVirtualSize(pan_netw.GetSize());
    
    # draw edges
    n=graph["node"];
    idl=graph["id2lab"];
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
        if len(lp) > 2:
            canvas.AddSpline(lp);
        else:
            canvas.AddLine(lp);
    for pn in graph["node"].values():
        if not (pn["label"] in pmetabs or pn["gr"]["type"]=="RECTANGLE"):
            continue;
        canvas.AddRectangle((pn["x"], pn["y"]), (pn["w"], pn["h"]));
    # draw nodes (rounded rectangles, rombs, etc.)
    # and their labels
    for pn in graph["node"].values():
        if pn["label"] not in pmetabs:
            if pn["gr"]["type"]=="ROUNDED_RECTANGLE":
                canvas.AddRectangle((pn["x"], pn["y"]), (pn["w"], pn["h"]), FillColor=pn["gr"]["fill"]);
            elif pn["gr"]["type"]=="DIAMOND":
                canvas.AddPolygon(dia2points(pn["x"], pn["y"], pn["w"], pn["h"]), FillColor=pn["gr"]["fill"]);
        else:
            if pn["label"] in netan["input"]:
                # inputs are in green
                canvas.AddRectangle((pn["x"], pn["y"]), (pn["w"], pn["h"]), LineColor="Green", FillColor="White");
            else:
                # regular nodes
                canvas.AddRectangle((pn["x"], pn["y"]), (pn["w"], pn["h"]), FillColor="White");
        canvas.AddText(pn["label"],
            (pn["x"]+0.5*pn["w"], pn["y"]+0.5*pn["h"]), Size=10, Position="cc");
    canvas.ZoomToBB();
    canvas.Draw(True);

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
        n["y"]=-float(n["gr"]["y"]);
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
            item["y"]=-float(item["y"])+dy;
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
    w*=0.5;
    h*=0.5;
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

def drawBars(ti_m):
    """drawBars(ti_m) -> None
    Draw rectangle bars in the metabolite nodes representing
    labeling propagation.
    """
    global bars;
    nds=graph["node"];
    nb_bars=data.shape[1];
    if not bars:
        # prepare list of tuples (x,y,w) of lower left corner of each rectangle
        bars=[None]*nb_bars;
        for m in pmetabs:
            if m not in dind:
                continue; # no data for this metabolite, it must be an input
            nd=nds[m];
            py=nd["y"];
            cols=dind[m];
            n=len(cols);
            w=nd["w"]/n;
            pxs=[round(nd["x"]+i*w) for i in xrange(n)];
            pxs.append(nd["x"]+nd["w"]);
            #print("m=", m, "pxs=", pxs);
            for (i,ic) in enumerate(cols):
                x,y,w=(pxs[i],py,pxs[i+1]-pxs[i]);
                color=wx.NamedColour(settings["colour_bar_bold"]);
                r,g,b=color.Get();
                color.Set(r,g,b, alpha=128);
                r=MovingRectangle((x,y), (w,0), FillColor=color, LineWidth=0);
                r.PutInForeground();
                bars[ic]=r;
                canvas.AddObject(r);
    #print("bars=", bars);
    #print("n=", [(nds[l]["x"], nds[l]["y"], nds[l]["w"], l) for l in pmetabs]);
    #sys.exit(1);
    # find the max row index ir in data s.t. ti(ir)<=tm
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
    # set the heights of all bars
    # the bars ar proportional to the faction of labeled component
    #print("ti_m=", ti_m, "ir=", ir, "f=", f);
    for m in pmetabs:
        if m not in dind:
            continue; # no data for this metabolite, it must be an output
        nd=nds[m];
        cols=dind[m];
        H=nd["h"];
        for ic in cols:
            v=data[ir,ic]*(1.-f)+data[ir+1,ic]*f;
            r=bars[ic];
            x,y,w,h=r.GetXYWH();
            r.SetShape((x,y), (w,round(H*v)));
            #print("ic=", ic, "r=", r);
            #r.Show(True);
    canvas.Draw(True);

def lab2ip(nb, lab):
    """lab2i(nb, nm) -> (i,Page) or (None,None)
    get page of a notebook nb by its label lab
    """
    nb_pages=nb.GetPageCount();
    for i in xrange(nb_pages):
        if lab == nb.GetPageText(i):
            return((i,nb.GetPage(i)));
    return((None,None));

def sett2wx(win):
    "create input wx widgets for settings in the window win"
    global sett_inp;
    gs=wx.FlexGridSizer(len(settings), 2, 2, 2)  # rows, cols, vgap, hgap;
    panel=wx.Panel(win, wx.ID_ANY);
    # create widgets and add them to the sizer
    gs.AddSpacer(10);
    gs.AddSpacer(10);
    sett_inp={};
    for (k,v) in settings.iteritems():
        #print("set k=", k);
        gs.Add(wx.StaticText(panel, wx.ID_ANY, k), 0, wx.EXPAND, 5);
        inp=wx.TextCtrl(panel, wx.ID_ANY, str(v), size=(125, -1));
        inp.Bind(wx.EVT_TEXT, OnSettings);
        inp.set_field=k;
        gs.Add(inp, 0, wx.EXPAND, 5);
        sett_inp[k]=inp;
    gs.AddSpacer(10);
    gs.AddSpacer(10);
    panel.SetSizer(gs);
    panel.Fit();
    win.SetVirtualSize(panel.GetSize());
    #panel.Center(wx.HORIZONTAL);

def show_tfactor(f):
    "set slider of speed, its label and setting input"
    settings["tw2tm"]=f;
    speed.SetValue((math.log(f)/math.log(4))*50+50);
    lf="%.1f"%f;
    labspeed.SetLabel("" if lf=="1.0" else "x"+lf if f > 1. else "/%.1f"%(1./f));
    if sett_inp["tw2tm"].GetValue() != str(f):
        sett_inp["tw2tm"].SetValue(str(f));

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
    wd=os.path.dirname(os.path.abspath(fftbl));
    os.chdir(wd);
    fftbl=os.path.basename(fftbl);
else:
    wd=os.path.abspath(os.getcwd());

## create GUI
app=wx.App();
code=tls.wxlay2py(tls.kvh2tlist(os.path.join(dirx, me[:-3]+"_lay.kvh")));
flname=sys.argv[0][:-3]+"_lay.py";
fp=open(flname, "w");
fp.write(code);
fp.close();
execfile(flname);
# add setting controls
sett2wx(sw_sett);

if fftbl:
    get_proj(os.path.join(os.getcwd(), fftbl));

app.MainLoop();

