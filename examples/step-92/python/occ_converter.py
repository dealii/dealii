#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
occ_converter  --- version 1.0, October 2024 ---

class  "occ_converter" for converting step-files into
 a) serialized (pickle-) objects
 b) gmsh .geo-files

(C) 2023-2024, Stephan Voss, stvoss@gmx.de, Neunkirchen am Brand, Germany

 * ---------------------------------------------------------------------
 *
 * This is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found at
 * https://www.gnu.org/licenses/old-licenses/lgpl-2.1
 * https://www.gnu.org/licenses/licenses/lgpl-3.0
 *
 * ---------------------------------------------------------------------


Actually, an "occ_converter" object is needed for each component of a model and
a step file may only represent a single component.

NOTE:
    This is actually a workaround and shall be further developed in the future such that step-files containing 
    complete assemblies and component-IDs should be handled at once 
    instead of having to create separate step-files for each model component
    (as it is in this example and the used classes in it's actual version in October 2024) 

"""
import gmsh     # used for importing OCC models (i.e. from step-files)
import math     # ...self-explaining ;-)
import sys      # used for some system- and file-operations
import pickle   # used for serialization, import and export as pickle-objects

"""
add an offset to entity numbers
"""
def ent_add_offset(_in_ent,_offset):
    return _in_ent + (_offset if _in_ent>=0 else -_offset)

# gmsh .geo entity names according to occ_converter's level number:
g_level_titles = ['Point',                  # level 0 (type 'Point')
                  'Line',
                  'Curve Loop',             # ...
                  'Plane Surface',          # ...
                  'Surface Loop',
                  'Volume',                 # level 5 (type 'Volume')
                  'Physical Volume',
                  'Physical Surface']

"""
The occ_converter class
"""
class occ_converter:
    
    Name = None         # This model's clear name (as string)
    
    pointdict = {}      # collection point entities
    linedict = {}       # collection line entities
    curvedict = {}      # collection curve entities
    surfacedict = {}    # collection surface entities
    volumedict = {}     # collection colume entities
    
    physical = None     # collection physicals
    
    entities = []
    
    point_entity_range = []
    line_entity_range = []
    curve_entity_range = []
    surface_entity_range = []
    
    geo_export_scale_factor = 1.0

    # object initializer
    def __init__(self,_geo_export_scale_factor = 1.0,_name=None):
        self.Name = _name
        self.clear()
        self.geo_export_scale_factor = _geo_export_scale_factor
    
    
    # reset all collections, dicts, ...
    def clear(self):
        
        self.Name = None
        self.entities.clear()
        
        self.pointdict.clear()
        self.linedict.clear()
        self.curvedict.clear()
        self.surfacedict.clear()
        self.volumedict.clear()
        self.physical = None

        self.point_entity_range = [None,None]
        self.line_entity_range = [None,None]
        self.curve_entity_range = [None,None]
        self.surface_entity_range = [None,None]
        
    
    # get the range of all used entity ID numbers (per type / level)    
    def calc_entities_min_max(self):
        # Points (level 0)
        for e,_ in self.pointdict.items():
            if self.point_entity_range[0]==None or self.point_entity_range[0]>e:
                self.point_entity_range[0] = e
            if self.point_entity_range[1]==None or self.point_entity_range[1]<e:
                self.point_entity_range[1] = e
        
        # Lines (level 1)
        for e,_ in self.linedict.items():
            if self.line_entity_range[0]==None or self.line_entity_range[0]>e:
                self.line_entity_range[0] = e
            if self.line_entity_range[1]==None or self.line_entity_range[1]<e:
                self.line_entity_range[1] = e
        
        # Curves (level 2)
        for e,_ in self.curvedict.items():
            if self.curve_entity_range[0]==None or self.curve_entity_range[0]>e:
                self.curve_entity_range[0] = e
            if self.curve_entity_range[1]==None or self.curve_entity_range[1]<e:
                self.curve_entity_range[1] = e
        
        # Surface Loops (level 3)
        for e,_ in self.surfacedict.items():
            if self.surface_entity_range[0]==None or self.surface_entity_range[0]>e:
                self.surface_entity_range[0] = e
            if self.surface_entity_range[1]==None or self.surface_entity_range[1]<e:
                self.surface_entity_range[1] = e
                
                
        return [self.point_entity_range[-1],self.line_entity_range[-1],self.curve_entity_range[-1],self.surface_entity_range[-1]]
        
        
    """
    sort (line / level 1) entities in a curve (level 2)
    
    in such way, that (signed) lines entities are listed in a correct sequence of points (entities) which are part of the curve 
    """    
    def sort_curve(self,_curve_entities):
        
        l_new_curve_entities = []
        l_curve_entities = [x[-1] for x in _curve_entities[:]]
        return l_curve_entities
    
        print('curve ',l_curve_entities)
        l_curve_entities_size = l_curve_entities.__len__()
        if l_curve_entities_size<=0:
            return l_new_curve_entities
        
        l_new_curve_entities.append(l_curve_entities[0])
        l_tmp_line = self.linedict[int(math.fabs(l_new_curve_entities[-1]))][::1] if l_new_curve_entities[-1]>=0 else self.linedict[int(math.fabs(l_new_curve_entities[-1]))][::-1]
        print('adding ',l_new_curve_entities[-1],' --> ', l_tmp_line)
        
        del l_curve_entities[0]
        
        l_curve_entities_size -= 1
        l_oriented_line_entity = l_new_curve_entities[-1]
        l_line_entity = l_oriented_line_entity if l_oriented_line_entity>=0 else -l_oriented_line_entity
        l_contact_entity = self.linedict[l_line_entity][-1] if l_oriented_line_entity>=0 else self.linedict[l_line_entity][0]
        
        while l_curve_entities.__len__()>0:
            for l_ix in range(l_curve_entities_size):
                l_tmp_oriented_line_entity = l_curve_entities[l_ix]
                l_tmp_line_entity = l_tmp_oriented_line_entity if l_tmp_oriented_line_entity>=0 else -l_tmp_oriented_line_entity
                l_cond_pos = l_contact_entity == self.linedict[l_tmp_line_entity][0] 
                l_cond_neg = l_contact_entity == self.linedict[l_tmp_line_entity][-1] 
                if l_cond_pos or l_cond_neg:
                    l_new_curve_entities.append(l_tmp_line_entity if l_cond_pos else -l_tmp_line_entity)
                    l_tmp_line = self.linedict[int(math.fabs(l_new_curve_entities[-1]))] if l_cond_pos else self.linedict[int(math.fabs(l_new_curve_entities[-1]))][::-1]
                    print('adding ',l_new_curve_entities[-1],' --> ', l_tmp_line)
                    
                    del l_curve_entities[l_ix]
                    
                    l_curve_entities_size -= 1
                    l_oriented_line_entity = l_new_curve_entities[-1]
                    l_line_entity = l_oriented_line_entity if l_oriented_line_entity>=0 else -l_oriented_line_entity
                    l_contact_entity = self.linedict[l_line_entity][-1] if l_oriented_line_entity>=0 else self.linedict[l_line_entity][0]
                    break
                    
        return l_new_curve_entities    
            
    
    """
    import a step-file utilizing gmsh
    and collect its entities (Points, Lines, Curves, ....) in this object
    """
    def occ_import(self,_infilename):
        # 1) use gmsh's OCC interface for importing the step-file 
        self.clear()
        
        self.Name = _infilename
        
        gmsh.initialize(interruptible=False)
        gmsh.clear()
        #gmsh.option.setNumber("General.NumThreads", 32)
        #gmsh.option.setNumber("Geometry.OCCParallel", 1)
        #print("nt:", gmsh.option.getNumber("General.NumThreads"))
        gmsh.model.add(_infilename)
        
        # Let's build the same model as in `t5.py', but using constructive solid
        # geometry.
        
        # We can log all messages for further processing with:
        gmsh.logger.start()
        
        # need to synchronize before looking up the point
        
        self.occ_entities = gmsh.model.occ.importShapes(_infilename)
        gmsh.model.occ.synchronize()
        
        # 2) transform the gmsh.occ model into collections of entities in this object
        self.entities = []
        for dim in range(4):
            self.entities.append( gmsh.model.occ.getEntities(dim) )
        
        no_parametrization=[]
        
        # Points
        for d,e in self.entities[0]:
            c = gmsh.model.getValue(0, e, no_parametrization)
            self.pointdict[e] = c
        
        # Lines
        for d,e in self.entities[1]:
            c = gmsh.model.get_boundary([[d,e]], False,True,False)
            self.linedict[e] = [int(c[x][-1]) for x in range(c.__len__())]
            ps = ['{0:.0f}'.format(x[1]) for x in c]

        # Curves        
        for d,e in self.entities[2]:
            c = gmsh.model.get_boundary([[d,e]], False,True,False)
            sc = self.sort_curve(c)
            self.curvedict[e] = [int(x) for x in sc]
        
        # Surface Loops
        for d,e in self.entities[3]:
            c = gmsh.model.get_boundary([[d,e]], False,True,False)
            ce = [int(x[1]) for x in c]
            self.surfacedict[e] = ce
            self.volumedict[e] =  [ e ]
        
        self.calc_entities_min_max()
        
        gmsh.finalize()
        

    """
    add offsets to this objects entity ID numbers
    
    This may be necessary to merge several "occ_converter" objects, where entity ID numbers may overlap.
    
    _offsets[0]    : offset to Point elements (level 0)
    _offsets[1]    : offset to Line elements (level 1)
    _offsets[2]    : offset to Curve elements (level 2)
    _offsets[3]    : offset to Surface Loop elements (level 3)
    """
    def offset_all_entities(self,_offsets):
        l_pointdict = self.pointdict.copy()
        l_linedict = self.linedict.copy()
        l_curvedict = self.curvedict.copy()
        l_surfacedict = self.surfacedict.copy()
        l_physical = self.physical

        self.clear()
        
        for e,c in l_pointdict.items():
            self.pointdict[ent_add_offset(e,_offsets[0])] = c
        
        for e,c in l_linedict.items():
            self.linedict[ent_add_offset(e,_offsets[1])] = [ent_add_offset(x,_offsets[0]) for x in c]
        
        for e,c in l_curvedict.items():
            self.curvedict[ent_add_offset(e,_offsets[2])] = [ent_add_offset(x,_offsets[1]) for x in c]
        
        for e,c in l_surfacedict.items():
            self.surfacedict[ent_add_offset(e,_offsets[3])] = [ent_add_offset(x,_offsets[2]) for x in c]
            
        self.calc_entities_min_max()
        
        self.physical = l_physical
        
    """
    add offsets to this objects entity ID numbers by taking the maximum entity ID numbers from _other object
    """
    def offset_all_entities_by_other(self,_other):
        l_ent_offsets = [ _other.point_entity_range[-1], _other.line_entity_range[-1], _other.curve_entity_range[-1], _other.surface_entity_range[-1] ]
        self.offset_all_entities(l_ent_offsets)
        
    
    """
    geometrically scale this object (by multiplying a scale "_factor" to points coordinates
    """
    def scale_all_points(self,_factor):
        
        for e,c in self.pointdict.items():
            self.pointdict[e] = [x*_factor for x in c]
        
    """
    find the closes point to "_p"
    """
    def find_closest_point(self,_p,_tol=1e-3):
        rp = None
        min_dist = _tol
        
        for e,c in self.pointdict.items():
            
            l_d = [_p[ix]-c[ix] for ix in range(3) ]
            l_d = [ x*x for x in l_d ]
            
            l_dist = math.sqrt(sum(l_d))
            if l_dist<min_dist:
                rp = [e,c]
                min_dist = l_dist
            
        return rp
     
    """
    find a matching line between points "_p1" and "_p2"
    consider even the orientation of the line, which for positive orientation shall start in "_p1" and end in "_p2"
    """
    def find_matching_line(self,_p1,_p2):
        for e,c in self.linedict.items():
            
            if _p1==c[0] and _p2==c[1]:
                return e
            elif _p1==c[1] and _p2==c[0]:
                return -e
            
        return None
     
    """
    return all entities, where Point _p is contained
    """ 
    def find_point_in_higher_order(self,_p):
        # dict for return-value:
        rv = {1:[],       # level 1: Lines
              2:[],       # level 2: Curves 
              3:[]}       # level 3: Surface Loops
        
        for e,c in self.linedict.items():
            a_c = [(x if x>=0 else -x) for x in c]
            if _p in a_c:
                rv[1].append(e)
        
        for e in rv[1]:
            for c_e,c_c in self.curvedict.items():
                a_c = [(x if x>=0 else -x) for x in c_c]
                if e in a_c and c_e not in rv[2]:
                    rv[2].append(c_e if c_e>=0 else -c_e)
        
        for e in rv[2]:
            for s_e,s_c in self.surfacedict.items():
                a_c = [(x if x>=0 else -x) for x in s_c]
                if e in a_c and s_e not in rv[3]:
                    rv[3].append(s_e if s_e>=0 else -s_e)
                    
        return rv
         
    
    """
    return all entities (level 0,1,2,3), where this object is in contact with _other object
    
    return is a dict of zipped entity lists (this objects entities, other object entities)
    """
    def find_contact_surfaces(self,_other):
        
        found = []
        
        for other_p_ent,other_p_coord in _other.pointdict.items():
            tmp_found = self.find_closest_point(other_p_coord)
            if tmp_found!=None:
                found.append([tmp_found[0],tmp_found[1],other_p_ent,other_p_coord])
                
        
        l_match = {0:[], 1:[], 2:[],3:[]}

        match_from_this = None
        match_from_other = None
        
        for p_ix in range(found.__len__()):
            p = found[p_ix]
            if p_ix==0:
                match_from_this =  self.find_point_in_higher_order(p[0])
                match_from_other =  _other.find_point_in_higher_order(p[2])
                
                match_from_this[0] = [ p[0] ]
                match_from_other[0] = [ p[2] ]
            else:
                match_from_this[0].append( p[0] )
                match_from_other[0].append( p[2] )

                tmp_match_from_this =  self.find_point_in_higher_order(p[0])
                tmp_match_from_other =  _other.find_point_in_higher_order(p[2])
                
                match_from_this[1].extend(tmp_match_from_this[1])
                match_from_other[1].extend(tmp_match_from_other[1])
                
                match_from_this[2] = list(set(match_from_this[2]).intersection(tmp_match_from_this[2]))
                match_from_other[2] = list(set(match_from_other[2]).intersection(tmp_match_from_other[2]))
                
        if match_from_this!=None:
            #match_from_this[1] = list(set(match_from_this[1]))
            #match_from_other[1] = list(set(match_from_other[1]))
        
            l_match[0] = list(zip(match_from_this[0],match_from_other[0]))
            l_match[1] = list(zip(match_from_this[1],match_from_other[1]))
            l_match[2] = list(zip(match_from_this[2],match_from_other[2]))
            l_match[3] = list(zip(match_from_this[3],match_from_other[3]))
    
        return l_match


    """
    In a list "lst" of entities, replace all occurrences of entity ID "old" by entity ID "new" 
    """
    @staticmethod    
    def list_replace(lst, old=1, new=10):
        """replace list elements (inplace)"""
        i = -1
        try:
            while True:
                i = lst.index(old, i + 1)
                lst[i] = new
        except ValueError:
            pass
        
    
    """
    In this objects Lines, replace all occurrences of the points listed 
    in a first index in the _pointzip zip by the second index in the _pointzip zip
    """
    def replace_points_in_lines(self,_pointzip):
        for point_pair in _pointzip:
            self.list_replace(self.linedict,point_pair[0],point_pair[1])
    
    
    """
    get_ext_surface_and_interface_curves():
    
    return two lists of entities:
    1) r_surface : list of surface entities without contact to any other surfaces
    2) r_interface: list of "surface" entities in contact with other surfaces
    
    The difference is simply:
    - interface curves will occur more than 1 time in the surfacedict
    - real surfaces only occur 1 single time in the surfacedict
    """
    def get_ext_surface_and_interface_curves(self):
        l_count_dict={}
        for s_e,s_c in self.surfacedict.items():
            for l_c in s_c:
                tmp_c = int(math.fabs(l_c))
                if l_count_dict.get(tmp_c,None)==None:
                    l_count_dict[tmp_c] = 1
                else:
                    l_count_dict[tmp_c] += 1
                    
        r_surface = []
        r_interfaces = []
        
        for tmp_c,counter in l_count_dict.items():
            if  counter<2:
                r_surface.append(tmp_c)
            else:
                r_interfaces.append(tmp_c)
                
                
        return r_surface, r_interfaces
        
    
    """
    merge all entities from _other object into this one
    
    This routine is keeping clean the entities lists of this and _other objects
    
    "merge()" function simply shuffles the lists from this and _other object into joined lists,
    but isn't doing any cleaning-up.
     
    """
    def merge_entities(self,_other):
        """
        - assign new surface entity IDs
        - for evesy surface entity ID in _other:
            - replace the touching contour ID in _other by _this
            - out of touching contour ID in _other:
              - replace every replace every line ID in _other by _this
        """
        
        #_other.offset_all_entities_by_other(self)
        #_other.geo_export()
        
        
        l_match_dict = self.find_contact_surfaces(_other)
        l_match_dict_unzipped = {}
        l_match_dicts_by_other = {}
        
        for key,val in l_match_dict.items():
            if val.__len__()>0:
                a,b = list(zip(*val))
                l_match_dict_unzipped[key] = a,b

                l_match_dicts_by_other[key] = dict(zip(b,a))
            
            else:
                l_match_dict_unzipped[key] = [],[]
            
                l_match_dicts_by_other[key] = {}
        
        l_matching_lines_other = []
        
        for matching_surface in l_match_dict[3]:
            """
            for every surface containing matching curves:
            - replace the matching curve in _other by the one from _this
            - delete the matching curve entities in _other.curvedict
            """

            for l_ix in range(l_match_dict_unzipped[2][-1].__len__()):
                l_this_ent  = int(math.fabs(l_match_dict_unzipped[2][0][l_ix]))
                l_other_ent = int(math.fabs(l_match_dict_unzipped[2][-1][l_ix]))
                #print('  replace plane ',l_other_ent,'by',l_this_ent,'in plane surface',matching_surface[-1])
                occ_converter.list_replace(_other.surfacedict[matching_surface[-1]],l_other_ent,l_this_ent)
                occ_converter.list_replace(_other.surfacedict[matching_surface[-1]],-l_other_ent,-l_this_ent)
            
                l_matching_lines_other.extend([int(math.fabs(x)) for x in _other.curvedict[l_other_ent]])
                #print('  cutting plane surface',l_other_ent)
                del _other.curvedict[l_other_ent]
                
        l_matching_lines_other = list(set(l_matching_lines_other))
        
        for key,val in l_match_dicts_by_other[1].items():
            """
            for line containing point matches:
            - replace the matching points in _other by the one from _this
            """
            if _other.linedict.get(key,None)!=None:
                for l_other_ent in l_match_dict_unzipped[0][-1]:
                    l_this_ent = l_match_dicts_by_other[0][l_other_ent]
                    
                    occ_converter.list_replace(_other.linedict[key],l_other_ent,l_this_ent)
                    occ_converter.list_replace(_other.linedict[key],-l_other_ent,-l_this_ent)
                    

        """
        delete all matching points in _other
        """
        for matching_point in l_match_dict_unzipped[0][-1]:
            del _other.pointdict[matching_point]
        
        
        l_matching_line_dict = {}
        for l_matching_line_other in l_matching_lines_other:
            l_p = _other.linedict[l_matching_line_other]
            l_matching_line_dict[l_matching_line_other] = self.find_matching_line(l_p[0],l_p[1])
             
        #print(l_matching_line_dict)
            
            
        
        if True:
            for l_other_ent in l_matching_lines_other:
                l_this_ent = l_matching_line_dict[l_other_ent]
                
                #print('curve',key,'replace line ',l_other_ent, 'by', l_this_ent )
                
                
            for key,val in _other.curvedict.items():
                """
                replace all matching lines in all curves in _other:
                """
                for l_other_ent in l_matching_lines_other:
                    #l_this_ent = l_match_dicts_by_other[1][l_other_ent]
                    l_this_ent = l_matching_line_dict[l_other_ent]
                    occ_converter.list_replace(_other.curvedict[key],l_other_ent,l_this_ent)
                    occ_converter.list_replace(_other.curvedict[key],-l_other_ent,-l_this_ent)
                
                    #print('curve',key,'replace line ',l_other_ent, 'by', l_this_ent )
                
            for l_matching_line_other in l_matching_lines_other:
                """
                for every matching curve:
                - delete all its lines entities in _other.linedict
                """
                del _other.linedict[l_matching_line_other]
            

            
        #print(l_match_dicts_by_other)


    """
    "merge()" function
    
    simply shuffles the lists from this and _other object into joined lists,
    but isn't doing any cleaning-up.
    
    (use merge_entities instead)
    
    This function may be useful for finding real surfaces and contact interfaces.
    """
    def merge(self,_other):
        self.pointdict.update(_other.pointdict)
        self.linedict.update(_other.linedict)
        self.curvedict.update(_other.curvedict)
        self.surfacedict.update(_other.surfacedict)
        self.calc_entities_min_max()
        

    """
    delete a Curve from this object (and its surfaces)
    """
    def delete_curve(self,_curve):
        del self.curvedict[_curve]

        try:
            while True:
                i = self.surfacedict.index(_curve, 0)
                del self.surfacedict[i]
                
        except ValueError:
            pass
        

    """
    get_geo_lines()
    
    get lists of strings for the levels of entities (see g_level_titles)
    
    if "_physical_only" is True: only level 4 and higher are exported (lists up to level 3 remain empty),
    because only those levels can be assigned a physical ID
    
    """
    def get_geo_lines(self,_physical_only=False):
        
        points=[]           # level 0
        lines=[]
        curves=[]           # ...
        surfaces=[]         # ...
        surface_loops=[]    
        volumes=[]          # level 5
        physicals=[]
        
        if not _physical_only:
            for e,c in self.pointdict.items():
                points.append('Point({0:.0f}) = '.format(e)+'{'+' {0:f}, {1:f}, {2:f} '.format(c[0]*self.geo_export_scale_factor,
                                                                                       c[1]*self.geo_export_scale_factor,
                                                                                       c[2]*self.geo_export_scale_factor)+'};')
        
            for e,c in self.linedict.items():
                ps = ['{0:.0f}'.format(x) for x in c]
                lines.append('Line({0:.0f}) = '.format(e)+'{ '+', '.join(ps)+' };')
            
            for e,c in self.curvedict.items():
                ps = ['{0:.0f}'.format(x) for x in c]
                curves.append('Curve Loop({0:.0f}) = '.format(e)+'{ '+', '.join(ps)+' };')
                surfaces.append('Plane Surface({0:.0f}) = '.format(e)+'{ '+'{0:.0f}'.format(e)+' };')
        
        physicals_strings=[]

        for e,c in self.surfacedict.items():
            c = list(set(c))
            ps = ['{0:.0f}'.format(x if x>=0 else -x) for x in c]
            slstring = ', '.join(ps)
            vstring = '{0:.0f}'.format(e)
            if not _physical_only:
                surface_loops.append('Surface Loop({0:.0f}) = '.format(e)+'{ '+slstring+' };')
                volumes.append('Volume({0:.0f}) = '.format(e)+'{ '+vstring+' };')

            if self.physical != None:
                if self.physical[0]==2:
                    physicals_strings.append(slstring)
                elif self.physical[0]==3:
                    physicals_strings.append(vstring)

        if physicals_strings.__len__()>0 and self.physical != None:
            fstring = '"{0:s}", {1:.0f}'.format(self.physical[1],self.physical[2]) if self.physical[2]!=0 else '0'
            if self.physical[0]==2:
                physicals.append('Physical Surface('+fstring+') = { '+', '.join(physicals_strings) +' };')
            elif self.physical[0]==3:
                physicals.append('Physical Volume('+fstring+') = { '+', '.join(physicals_strings) +' };')

        return points,lines,curves,surfaces,surface_loops,volumes,physicals


    """
    export this model as gmsh .geo-file
    """
    def geo_export(self,_output_file=sys.stdout,_physical_only=False):
        
        l_ofile = sys.stdout if _output_file==None else _output_file
        
        groups = self.get_geo_lines(_physical_only = _physical_only)
        
        for group in groups:
            for line in group:
                print(line,file=l_ofile)


    """
    export this serialized (pickle-) object
    """ 
    def pickle_export(self,_output_file):
        
        pickle.dump(self.pointdict, _output_file)
        pickle.dump(self.linedict, _output_file)
        pickle.dump(self.curvedict, _output_file)
        pickle.dump(self.surfacedict, _output_file)

    """ 
    export serialized (pickle-) object
    """ 
    def pickle_import(self,_input_file):
        
        self.Name = _input_file.name
        self.pointdict = pickle.load(_input_file)
        self.linedict = pickle.load(_input_file)
        self.curvedict = pickle.load(_input_file)
        self.surfacedict = pickle.load(_input_file)
        self.calc_entities_min_max()
