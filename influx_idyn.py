#!/usr/bin/env python
"""influx_idyn.py produces dynamic visualisation of isotopomer labeling
propagation in a metabolic network

usage: influx_idyn.py [-h|--help] [--DEBUG] [network[.ftbl]]
"""
# 2009-03-24 sokol
# Copyright 2009, MetaSys/INRA, France

# This file content:
# -imports
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
import math;

import wx;
import wx.html;
import wx.grid;
import wx.lib.plot as wxp;

import webbrowser;
import xml.dom.minidom as xd;

#sys.path.append('/home/sokol/dev/python/pytools');
#sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
dirx=os.path.dirname(sys.argv[0]) or os.getcwd();
sys.path.append(dirx);
import tools_ssg as tls;
#from C13_ftbl import *;
import C13_ftbl as c13;

## config constants
# program name
me=os.path.basename(sys.argv[0]);
# border width in pixels for widgets
wd_border=3;
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
    "frame_period": 1000, # milliseconds between two updates of the picture
};
#pnodes={}; # node plot attributes
pmetabs=set(); # set of plotted metabolites
dom=None;
#xnodes=None;
netw_img=None; # static image of the network
graph={}; # dictionary describing metab graph (node, edge and so on)
ti_wcs=None; # wall clocktime at the model's start
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
    pan_netw.SetBackgroundStyle(wx.BG_STYLE_CUSTOM);
    pbuf=netw_img.ConvertToBitmap();
    dc = wx.BufferedPaintDC(pan_netw, pbuf);

def OnLinkClicked(evt):
    #print(str(dir(evt)));
    #for o in dir(evt.GetLinkInfo()):
    #    print(o+"="+str(eval("evt.GetLinkInfo()."+o)));
    webbrowser.open_new_tab(evt.GetLinkInfo().Href);

def OnTimer(evt):
    """plot the current state of labeled matabs
    if evt==None, set timers to zero.
    """
    global ti_wcs, netw_img;
    if ti_wcs is None:
        ti_m=data[0,0]; # starting value for modeled time
        ti_wcs=time.time();
    ti_wcc=time.time();
    ti_m=(ti_wcc-ti_wcs)*settings["tw2tm"];
    ftime=s2ftime(ti_wcc-ti_wcs);
    pbuf=netw_img.ConvertToBitmap();
    dc=wx.BufferedPaintDC(pan_netw, pbuf);
    #dc.DrawLabel(ftime, (ti_wcc-ti_wcs, ti_wcc-ti_wcs, 100, 100));
    #print("ftime=", ftime);
    drawBars(dc, ti_m);
    if ti_m > data[-1,0]:
        # we are at the end of modeled time
        mainframe.timer.Stop();
        netw_img=pbuf.ConvertToImage();

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
    global network, pnetw, cnames, data, pan_data, dom, netan, pmetabs,\
        netw_img, dind, graph;
    if not fn or fn[-5:] != ".ftbl":
        return;
    # parse ftbl to netan
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
    # fill the set of metabolite names
    pmetabs=set(re.split("[:#+]", s)[0] for s in cnames[1:]);
    pmetabs.update(netan["input"]);
    # store indeces of data columns for each metabolite
    dind={};
    for (ic,nm) in enumerate(cnames):
        m=re.split("[:#+]", nm)[0];
        dind[m]=dind.get(nm, []);
        dind[m].append(ic);
    
    ## create and fill the grid
    if "grid" in dir(pan_data):
        # destroy the previous grid
        pan_data.grid.Destroy();
    # create new one
    pan_data.grid=wx.grid.Grid(pan_data, wx.ID_ANY);
    g=pan_data.grid;
    g.CreateGrid(data.shape[0], data.shape[1]);
    # fill the grid
    for ic in xrange(data.shape[1]):
        g.SetColLabelValue(ic, cnames[ic])
        for ir in xrange(data.shape[0]):
            g.SetCellValue(ir, ic, str(data[ir,ic]));
            g.SetReadOnly(ir, ic, True);
    # autosize the grid
    g.AutoSize();
    #print("gz=", g.GetSize());
    pan_data.SetClientSize(g.GetSize());
    sw_data.SetVirtualSize(pan_data.GetSize());
    
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
    pbuf=wx.EmptyBitmap(w+2*mrg_x, h+2*mrg_y);
    # connect the buffer to DC
    dc=wx.BufferedPaintDC(pan_netw, pbuf);
    # draw to buffer via DC
    drawGraph(graph, dc);
    netw_img=pbuf.ConvertToImage();
    #print(type(pbuf));
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
    pan_netw.SetSize((w+2*mrg_x, h+2*mrg_y));
    sw_netw.SetVirtualSize(pan_netw.GetSize());
    #pan_netw.Center();
    dc.Clear();
    dc.SetLogicalOrigin(-mrg_x,-mrg_y);
    # set font size
    font=dc.GetFont();
    font.SetPointSize(int(settings["fontPointSize"]));
    dc.SetFont(font);
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
            lp=lp[0:1]+\
                [(float(i["x"]), float(i["y"])) for i in b]+\
                lp[-1:];
        dc.DrawSpline(lp);
    # draw strait rectangles
    #print("nodes (x,y,w,h)=", [(pn["x"], pn["y"], pn["w"], pn["h"]) for pn in pnodes if pn["shape"]=="RECTANGLE"]);
    dc.DrawRectangleList([(pn["x"], pn["y"], pn["w"], pn["h"])
        for pn in graph["node"].values() if pn["label"] in pmetabs or pn["gr"]["type"]=="RECTANGLE"]);
    # draw other nodes (rounded rectangles, rombs, etc.) (not animated)
    # and their labels
    for pn in graph["node"].values():
        if pn["label"] not in pmetabs:
            if pn["gr"]["type"]=="ROUNDED_RECTANGLE":
                dc.SetBrush(wx.Brush(pn["gr"]["fill"]));
                dc.DrawRoundedRectangle(pn["x"], pn["y"], pn["w"], pn["h"], settings["radius"]);
            elif pn["gr"]["type"]=="DIAMOND":
                dc.SetBrush(wx.Brush(pn["gr"]["fill"]));
                #dc.DrawRectangle(pn["x"], pn["y"], pn["w"], pn["h"]);
                dc.DrawPolygon(dia2points(pn["x"], pn["y"], pn["w"], pn["h"]));
        dc.DrawLabel(pn["label"],
            (pn["x"], pn["y"], pn["w"], pn["h"]), wx.ALIGN_CENTER);

def dom2gr(dom):
    """dom2gr(dom)-> dict
    extract nodes and edges from dom structure (xgmml)
    """
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
    res["box"]=(minx, miny, w, h);
    print("gr w,h=", w, h);
    # get graph edges
    res["edge"]=[el2dict(e) for e in dom.getElementsByTagName("edge")];
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

def drawBars(dc, ti_m):
    """drawBars(dc, ti_m) -> None
    Draw rectangle bars in the metabolite nodes representing
    labeling propagation.
    """
    global base_rect;
    mrg_x=settings["mrg_x"];
    mrg_y=settings["mrg_y"];
    dc.SetLogicalOrigin(-mrg_x,-mrg_y);
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
    print("ti_m=", ti_m, "ir=", ir, "f=", f);
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
    dc.SetBrush(wx.Brush(settings["colour_bar_bold"]));
    dc.DrawRectangleList(rect);

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
        print("wd='"+wd+"'");
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
    nb_pages=nb.GetPageCount();
    for i in xrange(nb_pages):
        if sw_netw.Id == nb.GetPage(i).Id:
            nb.SetSelection(i); # show directly network tab
app.MainLoop();

