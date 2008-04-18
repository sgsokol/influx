# - Parse .ftbl (flux defs for C13flux)
# - Analyse ftbl

# 2008-01-22 sokol: ftbl_parse(f)
# 2008-01-25 sokol: ftbl_netan(ftbl)
# 2008-01-30 sokol: sto_r_m is added to the result of ftbl_netan()
# 2008-03-05 sokol: ftbl_parse(): added labeled atoms transitions in "TRANS" key
# 2008-03-05 sokol: ftbl_netan(): added isotopomers balance equations
# 2008-03-06 sokol: ftbl_netan(): added fwd-rev flux balance matrix
# 2008-03-10 sokol: ftbl_netan(): added carbon length
# 2008-03-?? sokol: ftbl_netan(): added chemical formula (reactions)
# 2008-03-?? sokol: ftbl_netan(): added cumomer balance pattern (cumo_balance_pattern)
# 2008-03-19 sokol: enum_metab(): metabolite list enumerated along to chemical pathways
# 2008-03-20 sokol: ftbl_netan(): added cumomer systems (cumo_sys)
# 2008-03-21 sokol: ftbl_netan(): added carbon transitions (carbotrans)
# 2008-03-25 sokol: iw(): added "index in given weight" iterator
# 2008-03-26 sokol: ftbl_netan(): added Cmax, maximal carbon length
# 2008-03-26 sokol: cumo_path(): enumerate cumomers of a given weight along pathways
# 2008-03-26 sokol: ftbl_netan(): removed metab_paths ordering
# 2008-03-26 sokol: ftbl_netan(): removed cumo_w ordering
# 2008-03-26 sokol: ftbl_netan(): added set of all metabolites (metabs)
# 2008-04-14 sokol: ftbl_netan(): added set of free fluxes (flux_free[net|xch])
# 2008-04-14 sokol: ftbl_netan(): added set of constrained fluxes flux_constr[net|xch])
# 2008-04-16 sokol: ftbl_netan(): added set of measured fluxes (flux_measured)
# 2008-04-18 sokol: ftbl_netan(): changed structure of sto_m_r, sto_r_m, flux_m_r, cumo_balance_pattern and matrix A

from tools_ssg import *;
def ftbl_parse(f):
    """reads .ftbl file. The only input parameter f is a stream pointer
    with read permission. It's user care to open and close it outside of
    this module. This function parses the input and returns a dictionnary
    with items corresponding to sections in .ftbl. One section is added.
    'TRANS' correponds to carbon transitions.""";
    import re;
    ftbl={};    # main dictionary to be returned;
    
##    print f;
    reblank=re.compile("^[\t ]*$");
    recomm=re.compile("^[\t ]*//.*$");
    reading="data";
    for l in f:
##        print "raw l="+l;
        # strip out \r
        l=l.replace("\r", "");
##        print "-ctrl-r l="+l;
        # strip out double quots
        if l[0] == '"': l=l[1:];    # very begining
##        print "-fq l="+l;
        if l[-1] == '"': l=l[:-1];    # very end
##        print "-lq l="+l;
        l=l.replace("\t\"", "\t");    # at the field begining
##        print "-ffq l="+l;
        l=l.replace("\"\t", "\t");    # at the end
##        print "-lfq l="+l;
        # skip comments and emty rows;
        if recomm.match(l) or reblank.match(l): continue;
        
        # skip the comments at the end of the row;
        l=re.sub(r'^(.+)//.*$',r'\1',l).rstrip();
        
        # split in fields;
        flds=l.split("\t");
        
##        print "proceeding:"+l;
        if reading=="data" and len(flds)==1:
            # new section starts here;
            sec_name=flds[0];
            subsec_name="";
            # prepare storage
            ftbl[sec_name]=[];
            # prepare storage for carbon transitions
            if sec_name=="NETWORK":
##               print "TRANS prepared"
               ftbl["TRANS"]={};
            try: del stock;
            except NameError: pass;
            stock=ftbl[sec_name];
            skiptab=1;
            data_count=0;
            reading="col_names";
            continue;
        if len(flds)==2 and len(flds[0])==0:
            # read subsection name;
            subsec_name=flds[1];
            # prepare sub-storage
            if not ftbl[sec_name]:
                ftbl[sec_name]={};
            ftbl[sec_name][subsec_name]=[];
            try: del stock;
            except NameError: pass;
            stock=ftbl[sec_name][subsec_name];
            skiptab=2;
            data_count=0;
            reading="col_names";
            continue;
        if reading=="col_names" and len(flds)>2:
            # read column names;
            col_names=l[skiptab:].split("\t");
            reading="data";
##            print "col_names=", col_names;
            continue;
        if reading=="data":
            data=l[skiptab:].split("\t");
            dic={};
            if sec_name == "NETWORK" and len(data[0])==0:
                # we are at carbon transition line here (e.g. #ABC -> #AB +#C)
##                print "data_count="+str(data_count), \
##                    "\ndata="+str([l for l in enumerate(data)]), \
##                    "\nstock="+str(stock);
                fl_name=str(stock[data_count-1][col_names[0]]);
                for i in range(len(col_names)):
                    try:
                        dic[col_names[i]]=data[i].strip();
                    except IndexError:
                        pass;
                ftbl["TRANS"][fl_name]=dic;
                continue;
            for i in xrange(len(col_names)):
                try:
                    dic[col_names[i]]=data[i].strip();
                except IndexError:
                    pass;
            stock.append(dic);
            data_count+=1;
##            print "sec, subsec=", sec_name, subsec_name, data_count, dic;
##        print "len(flds)=", len(flds), flds, l, data;
##        print "keys", ftbl.keys(), (ftbl[sec_name].keys() \
##                if type(ftbl[sec_name])==type({}) else "");
    return ftbl;

def ftbl_netan(ftbl):
    # analyse ftbl dictionary to find
    # - network inputs (input)
    # - network outputs (output)
    # - substrates (subs)
    # - products (prods)
    # - metabolites (metabs)
    # - reactions (reacs)
    # - not reversible reactions (subset of reacs) (notrev)
    # all above items are in named sets
    # - stocheometric matrix (sto_r_m)
    # - stocheometric matrix (sto_m_r)
    # - fwd-rev flux matrix (flux_m_r)
    # - cumomer balances (cumo_m_r_m)
    # - carbon length (Clen)
    # - reaction formula (formula)
    # - metabolite network (metab_netw)
    # - carbon transitions (carbotrans)
    # - free fluxes (flux_free)
    # - constrained fluxes (flux_constr)
    # - measured fluxes (flux_measured)

    # init named sets
    netan={
        'input':set(),
        'output':set(),
        'subs':set(),
        'prods':set(),
        'metabs':set(),
        'reac':set(),
        'notrev':set(),
        'sto_r_m':{},
        'sto_m_r':{},
        'flux_m_r':{},
        'cumo_balance_pattern':{},
        'Clen':{},
        'formula':{},
        'metab_netw':{},
        'cumo_sys':{},
        'carbotrans':{},
        'Cmax':{},
        'flux_free':{},
        'flux_constr':{},
        'flux_measured':{},
    };
    res="";     # auxiliary shot cut to current result;

    # analyse networks
    netw=ftbl.get('NETWORK');
    row_to_del=[];
    for row in netw:
        reac=row['FLUX_NAME'];
##        print "reac="+reac
        e1=row.get('EDUCT_1');
        if not e1: raise 'EDUCT_1 must be defined in flux '+reac;
        e2=row.get('EDUCT_2');
        p1=row.get('PRODUCT_1');
        if not p1: raise 'PRODUCT_1 must be defined in flux '+reac;
        p2=row.get('PRODUCT_2');
        
        # local substrate (es), product (ps) and metabolites (ms) sets
        es=set((e1,));
        if e2:
            es.add(e2);
        ps=set((p1,));
        if p2:
            ps.add(p2);
        ms=es|ps;
        
        if len(es)+len(ps)==4:
            # need of fictious metabolite AB (A+B=C+D => A+B=AB; AB=C+D)
            # add two reactions: fusion (f) and split (s)
            reac_f=reac+'.f';
            reac_s=reac+'.s';
            ab='.'.join(es);
            row_to_del.append(row);
            
            # add two new reactions in network and in flux conditions
            # and skip current reaction
            
            # fusion part
            ii=len(netw);
            netw.append({});
            netw[ii]['FLUX_NAME']=reac_f;
            netw[ii]['EDUCT_1']=e1;
            netw[ii]['EDUCT_2']=e2;
            netw[ii]['PRODUCT_1']=ab;
            # carbon transitions
            trans=ftbl['TRANS'][reac];
            ftbl['TRANS'][reac_f]={};
            ftbl['TRANS'][reac_f]['EDUCT_1']=trans['EDUCT_1'];
            ftbl['TRANS'][reac_f]['EDUCT_2']=trans['EDUCT_2'];
            ftbl['TRANS'][reac_f]['PRODUCT_1']=trans['EDUCT_1']+\
                trans['EDUCT_2'][1:];
            # flux constraints
            if 'FLUXES' in ftbl:
                if 'XCH' in ftbl['FLUXES'] and reac in ftbl['FLUXES']['XCH']:
                    ftbl['FLUXES']['XCH'][reac_f]=dict(ftbl['FLUXES']['XCH'][reac]);
                    ftbl['FLUXES']['XCH'][reac_s]=dict(ftbl['FLUXES']['XCH'][reac]);
                    del(ftbl['FLUXES']['XCH'][reac]);
                if 'NET' in ftbl['FLUXES'] and reac in ftbl['FLUXES']['NET']:
                    ftbl['FLUXES']['NET'][reac_f]=dict(ftbl['FLUXES']['NET'][reac]);
                    ftbl['FLUXES']['NET'][reac_s]=dict(ftbl['FLUXES']['NET'][reac]);
                    del(ftbl['FLUXES']['NET'][reac]);
            # split part
            ii=len(netw);
            netw.append({});
            netw[ii]['FLUX_NAME']=reac_s;
            netw[ii]['EDUCT_1']=ab;
            netw[ii]['PRODUCT_1']=p1;
            netw[ii]['PRODUCT_2']=p2;
            # carbon transitions
            trans=ftbl['TRANS'][reac];
            ftbl['TRANS'][reac_s]={};
            ftbl['TRANS'][reac_s]['EDUCT_1']=trans['PRODUCT_1']+\
                trans['PRODUCT_2'][1:];
            ftbl['TRANS'][reac_s]['PRODUCT_1']=trans['PRODUCT_1'];
            ftbl['TRANS'][reac_s]['PRODUCT_2']=trans['PRODUCT_2'];
            del(ftbl['TRANS'][reac]);
            continue;

        # normal reaction A+B=C or C=A+B
        netan['reac'].add(reac);
        netan['subs'].update(es);
        netan['prods'].update(ps);
        netan['metabs'].update(ms);
##        aff('ms for '+reac, ms);#

        # create chemical formula
        netan['formula'][reac]={'left':es, 'right':ps, 'all':ms};

        # Carbon length and transitions
        trans=ftbl['TRANS'][reac];
        netan['carbotrans'][reac]={'left':[], 'right':[]};
        for (m,carb,lr) in [(e1,trans.get('EDUCT_1'),'left'),
                (e2,trans.get('EDUCT_2'),'left'),
                (p1,trans.get('PRODUCT_1'),'right'),
                (p2,trans.get('PRODUCT_2'),'right')]:
##                print "m="+str(m), "; carb="+str(carb);
            if not m or not carb:
                continue;

            # carbon transitions
            netan['carbotrans'][reac][lr].append((m,carb[1:])); # strip '#' character

            # carbon length
            if netan['Clen'].get(m, 0) and \
                    netan['Clen'][m] != len(carb)-1:
                raise 'Metabolite '+m+' has length '+\
                        str(netan['Clen'][m])+' but in reaction '+reac+\
                        ' it has legth '+str(len(carb)-1);
            netan['Clen'][m]=len(carb)-1; # don't count '#' character

        # stocheometric matrix in dictionnary form
        # sto_r_m[reac]['left'|'right']=list(metab) and
        # sto_m_r[metab]['left'|'right']=list(reac)
        # substrates are in 'left' list
        # and products are in the 'right' one.

        netan['sto_r_m'][reac]={'left':[e1], 'right':[p1]};
        if (e2):
            netan['sto_r_m'][reac]['left'].append(e2);
        if (p2):
            netan['sto_r_m'][reac]['right'].append(p2);
##        print "col keys", col.keys();
        
        for s in netan['sto_r_m'][reac]['left']:
##           print "sto_m_r s="+str(s);
            if not s in netan['sto_m_r']:
                netan['sto_m_r'][s]={'left':[], 'right':[]};
            netan['sto_m_r'][s]['left'].append(reac);
        for s in netan['sto_r_m'][reac]['right']:
##           print "sto_m_r s="+str(s);
            if not s in netan['sto_m_r']:
                netan['sto_m_r'][s]={'left':[], 'right':[]};
            netan['sto_m_r'][s]['right'].append(reac);
        # end netan['sto_m_r']
##        aff('sto_m_r'+str(reac), netan['sto_m_r']);
    # cleanup splited and fused reactions
    for row in row_to_del:
        netw.remove(row);

    # find input and output
    netan['input'].update(netan['subs']-netan['prods']);
    netan['output'].update(netan['prods']-netan['subs']);

    # analyse reaction reversibility
    # and free fluxes
    # and constrained fluxes (dictionary flux->value)
    netan['flux_free']={'net':set(), 'xch':set()};
    netan['flux_constr']={'net':{}, 'xch':{}};
    flx=ftbl.get('FLUXES');
    if flx:
        xch=flx.get('XCH');
        net=flx.get('NET');
##        aff("net", net);
    else:
        sys.stdout=sys.stderr;
        print sys.argv[0]+': netan[FLUXES][XCH] is not defined';
        return None;
    for reac in netan['reac']:
        # reac can be .s (splitted) or .f (fused)
        # separate the name and sf_state
        tmp=reac.split('.');
        if (len(tmp) == 2):
            (oldreac, sf_state) = tmp;
        elif (len(tmp) == 1):
            # reac does not change
            oldreac=reac;
            sf_state='';
        else:
            raise 'Badly formatted reaction name <'+reac+'>';
        # get xch condition for this reac
        cond=[row for row in xch if row['NAME']==oldreac];
        # get net condition for this reac
        ncond=([row for row in net if row['NAME']==oldreac]) if net else [];
        if cond:
            cond=cond[0];
        if ncond:
            ncond=ncond[0];
        if (not cond and  not ncond):
            continue;
##        aff("cond", cond);
##        aff("ncond", ncond);
        # not reversibles are those reaction having xch flux=0 or
        # input/output fluxes
        if (cond and cond['FCD'] == 'C' and cond['VALUE(F/C)'] == '0') or \
                (netan['formula'][reac]['left'] & netan['input']) or \
                (netan['formula'][reac]['right'] & netan['output']):
            # input/output reactions cannot be .s or .f
            netan['notrev'].add(reac);
        elif (cond and cond['FCD'] == 'F'):
            # free xch
            netan['flux_free']['xch'].add(reac);
        elif (cond and cond['FCD'] == 'C'):
            # constrained xch
            netan['flux_constr']['xch'][reac]=float(cond['VALUE(F/C)']);
        if (ncond and ncond['FCD'] == 'F'):
            # free net
            netan['flux_free']['net'].add(reac);
        elif (ncond and ncond['FCD'] == 'C'):
            # constr net
            netan['flux_constr']['net'][reac]=float(cond['VALUE(F/C)']);
##        aff("free net", netan['flux_free']['net']);
##        aff("free xch", netan['flux_free']['xch']);

    # measured fluxes
    for row in ftbl.get('FLUX_MEASUREMENTS',[]):
        netan['flux_measured'][row['FLUX_NAME']]={\
                'value': row['VALUE'], \
                'deviation': row['DEVIATION']};
    
    # discard empty entries
    for e in netan:
        try:
            netan[e].discard('');
            netan[e].discard(None);
        except AttributeError: pass;

    # fwd-rev flux balance
    del(res);
    res=netan['flux_m_r'];
    # res[metab]['in']=list(flux.fwd|flux.rev)
    # res[metab]['out']=list(flux.fwd|flux.rev)
    for metab,lr in netan['sto_m_r'].iteritems():
        # lr is dico with 'left' and 'right' entries
##        print "fwd metab="+str(metab);
        if not metab in res:
            res[metab]={'in':[], 'out':[]};
        for reac in lr['left']:
            # here metabolite is consumed in fwd reaction
            # and produced in the reverse one.
            res[metab]['out'].append(reac+".fwd");
            if reac not in netan['notrev']:
                res[metab]['in'].append(reac+".rev");
        for reac in lr['right']:
            # here metabolite is consumed in rev reaction
            # and produced in the forward one.
            res[metab]['in'].append(reac+".fwd");
            if reac not in netan['notrev']:
                res[metab]['out'].append(reac+".rev");

    # cumomer balances (generic metabolite form)
    del(res);
    res=netan['cumo_balance_pattern'];
    # res[row_metab]['out']=flux_list;
    # res[row_metab]['in']=list( (flux,list(metab)) );
    for metab,iot in netan['flux_m_r'].iteritems():
        # for each metab the pattern is the same: "in" fluxes equals "out"
        # in 'in' part, we have tuples (flux,list(metabs))
        # only one metab corresponds to linear part
        # two cumomers correspond to quadratique part
        # in 'out' part, we have just a list of fluxes list(flux)
        res[metab]={'in':[], 'out':[]};
##        print "iso metab="+str(metab);
        for flux in iot['out']:
            # get output terms (supposed to be multiplied by cumomer itself)
            res[metab]['out'].append(flux);
        for flux in iot['in']:
            # get 'left' or 'right' for incoming metabolites
            lr=('left' if flux[-4:]=='.fwd' else 'right');
            reac=flux[:-4];
            # find flux index
            i=[fl for (fl,lst) in res[metab]['in'] if fl==flux];
            if (!len(i)):
                i=[len(res[metab]['in'])];
                res[metab]['in'].append(flux,[]);
            res[metab]['in'][i[0]][1].extend(netan['sto_r_m'][reac][lr]);

    # metabolite network
    del(res);
    res=netan['metab_netw'];
    # left part or reaction points to right part
    for reac,parts in netan['formula'].iteritems():
        for metab in parts['left']:
            if not metab in res:
                res[metab]=set();
            res[metab].update(parts['right']);
    # order cumomers by weight. For a given weight, cumomers are sorted by
    # metabolites order.
    netan['Cmax']=netan['Clen'][max(netan['Clen'])];
    Cmax=netan['Cmax'];
    
    # plain list of all metabs
    all_metabs=netan['metabs'];
    
    # cumomers systems A*x=b, one by weight
    # weights are going from 1 to Cmax where Cmax is the maximal
    # carbon string length in all metabolites
    # A is a list of matrices (in weight order)
    # b is a list of right hand parts (still in weight order)
    # the dimensions of various weights in b are not the same
    # too short metabolites are dropped when going to higher weights.
    del(res);
    res=netan['cumo_sys'];
    res['A']=[];
    res['b']=[];
##    aff('all', all_metabs);#
    for w in xrange(1,Cmax+1):
        # wheight w
        res['A'].append({});
        res['b'].append({});
        # run through metabs
        for metab in all_metabs:
            #print 'A keys 1:'+str(res['A'][w-1].keys()) ##
            # is this metab is long enough for this weihgt?
##            print metab;#
            Clen=netan['Clen'][metab];
            if w > Clen:
                # no equation for too short metabolite at this weight
                continue;
            if metab in netan['input'] or metab in netan['output']:
                # no equation for input neither output
                continue;
            
            # run through all cumomers for this metab
            # and linearize the system
            for iw in cumo_iw(w,Clen):
                row=netan['cumo_balance_pattern'][metab];
                # out (diagonal) terms
                # res['A'][w-1][row_cumo][row_cumo]=flux_list;
                flux_lst=row['out'].keys();
                m_cumo=metab+':'+str(iw);
##                print 'm_cumo='+m_cumo ##
##                print 'A keys 1.5:'+str(res['A'][w-1].keys()) ##
                # diagonal term
                res['A'][w-1][m_cumo]={m_cumo:flux_lst};
##                print 'A keys 2:'+str(res['A'][w-1].keys()) ##
                # in (off diagonal) terms
                # res['A'][w-1][row_cumo][col_cumo]=flux_lst
                for contribs,flux2coef in row['in'].iteritems():
                    for contrib in contribs.split(';'):
                        cumo2lflux={};
                        for flux,coef in flux2coef.iteritems():
                            flux_str=('' if coef==1 else str(coef)+'*')+flux;
                            reac=flux[:-4];
                            # calculate cumomer index of contributors from
                            # carbon transition relationships
                            carbotrans=netan['carbotrans'][reac];
                            metab_lr,contrib_lr=('left','right') \
                                if any(metab==m_s[0] for m_s in carbotrans['left']) \
                                else ('right','left');
                            # run through all occurences of metab in carbotrans
                            for cur_metab,mcarbostr in carbotrans[metab_lr]:
                                if cur_metab != metab:
                                    continue;
                                # run through all occurences of contrib
                                for cur_contrib,ccarbostr in carbotrans[contrib_lr]:
                                    if cur_contrib != contrib:
                                        continue;
                                    iwc=src_ind(ccarbostr, mcarbostr, iw);
                                    c_cumo=contrib+':'+str(iwc);
                                    wcontrib=sumbit32(iwc);
                                    if wcontrib != w:
                                        if wcontrib>0:
                                            # cumomer weight conservation is broken
                                            # drop this term in rigth hand side
                                            if not res['b'][w-1].get(m_cumo):
                                                res['b'][w-1][m_cumo]={};
                                            if not res['b'][w-1][m_cumo].get(flux):
                                                res['b'][w-1][m_cumo][flux]=set();
                                            # here we assume that flux coef==1
                                            res['b'][w-1][m_cumo][flux].add(c_cumo);
                                            continue;
                                        else:
                                            # it will become just a linear term
                                            continue;
                                    if not cumo2lflux.get(c_cumo):
                                        cumo2lflux[c_cumo]=[flux_str];
                                    else:
                                        cumo2lflux[c_cumo].append(flux_str);
                        for cumo, fluxes in cumo2lflux.iteritems():
                            if contrib in netan['input']:
                                for flux in fluxes:
                                    if not res['b'][w-1].get(m_cumo):
                                        res['b'][w-1][m_cumo]={};
                                    if not res['b'][w-1][m_cumo].get(flux):
                                        res['b'][w-1][m_cumo][flux]=set();
                                    # here we assume that flux coef==1 (are we right?)
                                    res['b'][w-1][m_cumo][flux].add(cumo);
                                continue;
                            res['A'][w-1][m_cumo][cumo]='-('+\
                                '+'.join(fluxes)+')';
                            #print 'A keys 3:'+str(res['A'][w-1].keys()) ##
                # end contribs,flux2coef in row['in']
    return netan;

def enum_path(starts, netw, outs, visited=set()):
    """Enumerate metabilites along to reaction pathways.
    Algo: start from an input, follow chemical pathways till an output or
    already visited metabolite. Returns a list of metabolite pathways.
    Each pathways is an ordered list."""
    res=[];
    if (len(visited) == len(netw.keys())):
        # no more metabolites 
        return res;
    for start in starts:
        if (start in visited):
           continue;
        visited.add(start);
        if start in outs:
            return [[start]];
        paths=enum_path(netw[start]-visited, netw, outs, visited);
        if (paths):
            res.append([start]+paths[0]);
            if len(paths) > 1:
                res.extend(paths[1:]);
        else:
            res.append([start]);
    return res;
def cumo_path(starts, A, visited=set()):
    '''Enumerate cumomers along to reaction pathways.
    Algo: start from an input, follow chemical pathways till no more
    neighbours or till rest only visited metabolite.
    Returns a list of cumomer pathways.
    Each pathways is an ordered list.'''
    #if [s for s in starts if s.find('.') >= 0]:
    #    aff('s', starts);
    res=[];
    if (len(visited) == len(A)):
        # no more cumomers
        return res;
    for start in starts:
        if (start in visited or not start in A):
           continue;
        visited.add(start);
        # next_starts are cumos influenced by start
        next_starts=set(cumo for cumo in A if start in A[cumo]);
##        aff('next to '+start, next_starts)
        next_starts-=visited;
        if next_starts:
            paths=cumo_path(next_starts, A, visited);
##            aff('paths for '+start, paths)
##            aff('visited', visited)
        else:
            path=[];
            res.append([start]);
            continue;
        if (paths):
            res.append([start]+paths[0]);
            if len(paths) > 1:
                res.extend(paths[1:]);
        else:
            res.append([start]);
##    aff('r', res)
    return res;
    
def src_ind(substrate, product, iprod):
    '''For given substrate and product carbon strings (e.g. "abc", "ab")
    calculate substrate index corresponding to product index'''
    movbit=1;
    isubstr=0;
    substrate=substrate[::-1];
    product=product[::-1];
    for nb in xrange(len(product)):
        if (movbit & iprod):
            # algo: if the current bit is set in product serach for its origin
            # and set the corresponding bit in substrate
            try:
                isubstr|=1<<substrate.find(product[nb]);
            except ValueError:
                pass;
        movbit<<=1;
    return isubstr;
def cumo_iw(w,nlen):
    '''iterator for a given cumomer wight w in the carbon length nlen'''
#    print ("cumo_iw: w=%d,nlen=%d\n" % (w,nlen))
    if w < 0 or w > nlen:
       raise "CumoWeightOutOfRange";
    if nlen < 0:
       raise "CumoLengthOutOfRange";
    if w == 1:
        movbit=1;
        for i in xrange(nlen):
            yield movbit;
            movbit<<=1;
    else:
        movbit=1<<(w-1);
        for i in xrange(nlen-w+1):
            for subi in cumo_iw(w-1,w+i-1):
                yield movbit+subi;
            movbit<<=1;
