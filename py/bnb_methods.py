import os 
import pandas as pd
import numpy as np
from effect_calculator import EffectCalculator
import pybnb
import math
import sys
import datetime
import graph_utils
import networkx as nx


# call this like: mpiexec -n 50 python bnb_methods.py PaperBnB 2

###
# Structure of the library:
# - Simple is the main BnB class, used to implement the other BnB approaches
# - PaperBnB is the reimplementation of the BnB approach proposed in "Discovering Reliable Causal Rules"
# - OurBnB is the version of BnB with the improved bound but IT WORKS ONLY ON SINGLE CORE since all the dictionaries (such as previous_layer_a_dict, current_layer_a_dict) must be sincronized between all the workers, but this is not implemented
# - ReliablePaperBnB is our version of BnB that outputs reliable rules with correction for multiple hypothesis testing
###

class Simple(pybnb.Problem):
    
    def __init__(self, data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, include_unmutated_treats = False):
        self.data_path = data_path
        self.cfnds = cfnds
        self.treats = treats
        self.treats_vals = treats_vals
        self.outcome = outcome
        self.outcome_vals = outcome_vals
        Simple.total_treats = total_treats
        Simple.max_rule_len = max_rule_len
        self.include_unmutated_treats = include_unmutated_treats
            
    def __repr__(self) -> str:
        s = EffectCalculator.rule_tostring(self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals) + "\n"
        return s

    def __eq__(self, __o: object) -> bool:
        res = True
        if not isinstance(__o, Simple) and not isinstance(__o, pybnb.Node):
            return False
        if isinstance(__o, Simple):
            print("YO")
            res = res and __o.treats == self.treats
            res = res and __o.treats_vals == self.treats_vals
        if isinstance(__o, pybnb.Node):
            print("YO1")
            res = res and __o.state[2] == self.treats
            res = res and __o.state[3] == self.treats_vals
        return res

    #
    # required methods
    #
    def sense(self):
        # return pybnb.minimize
        return pybnb.maximize
    
    #
    # optional methods
    #
    def notify_solve_begins(self,
                            comm,
                            worker_comm,
                            convergence_checker):
        pass

    # def notify_new_best_node(self,
    #                          node,
    #                          current):
    #     pass

    def notify_solve_finished(self,
                              comm,
                              worker_comm,
                              results):
        pass

    def notify_new_best_node(self,
                             node,
                             current):
        # Simple.CURR_BEST = max(Simple.CURR_BEST, EffectCalculator.rule_tostring(self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals))
        # print("Found new best node",Simple.CURR_BEST)
        # if current and node ==self:
        # # if current:
        #     print("A",self,"B",node.state,"C",current)
        pass


class PaperBnB(Simple):
    
    def __init__(self, data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, per_strata_data = [], include_unmutated_treats = False):
        super().__init__(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, include_unmutated_treats = include_unmutated_treats)
        ### PAY ATTENTION: at the init phase and on bound() funct, self.per_strata_data contains the strata infos of the parent
        # and after the call objective() it contains the strata infos of the current node (to be passed to the children)
        self.per_strata_data = per_strata_data
            
    def __repr__(self) -> str:
        s = super().__repr__()
        s += "len(self.per_strata_data):"+str(len(self.per_strata_data))
        return s

    #
    # required methods
    #
    
    def objective(self):
        self.df = pd.read_csv(self.data_path)
        if len(self.treats) == 0:
            return -math.inf
        
        # must keep the Laplace correction to work with the tight optimistic estimator
        # after the call self.per_strata_data contains the value of this node's strata and must be passed to child for optimistic bound estimation
        effect, self.per_strata_data = EffectCalculator.calculate_effect(self.df, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, laplace_correction = True, return_per_strata_data = True)
        return effect
    
    def bound(self):
        if len(self.treats) < 2:  # I don't have a value set up for per_strata_data
            return math.inf 
        optimist_effect = 0
        tot_cnt = 0
        # when I call bound() self.per_strata_data contains the value of each father's strata
        for pair in self.per_strata_data:
            a,b,c,d,n = pair
            n1 = a + c
            strata_best_eff = - math.inf
            for aa in range(a+1):
                curr_eff = (aa+1)/(aa + 2) - (n1-aa+1)/(n-aa+2) # this implicitely assumes Laplace correction
                if curr_eff > strata_best_eff:
                    strata_best_eff = curr_eff
            optimist_effect += strata_best_eff*n
            tot_cnt += n
        optimist_effect = optimist_effect/tot_cnt

        return optimist_effect
    
    def save_state(self, node):
        node.state = (self.data_path, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, self.per_strata_data)

    def load_state(self, node):
        (self.data_path, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, self.per_strata_data) = node.state

    def branch(self):
        if len(self.treats) + max(1,len(self.cfnds)) < PaperBnB.max_rule_len: # with the max I consider the emptyset as a set of dim. 1 (since we correct on all single Zs and the empty set by construction)
            for v in PaperBnB.total_treats:
                if v not in self.treats:
                    v_vals_list = [1]
                    if self.include_unmutated_treats:
                        v_vals_list.append(0)
                    for v_val in v_vals_list:
                        child = pybnb.Node()
                        if self.treats == []:
                            child_treats = [v]
                            child_treats_vals = [[v_val]]
                        else:
                            child_treats = self.treats + [v]
                            child_treats_vals = self.treats_vals + [[v_val]]
                        child.state = (self.data_path, self.cfnds, child_treats, child_treats_vals, self.outcome, self.outcome_vals, self.per_strata_data)
                        yield child
                
    
class OurBnB(Simple):
    previous_layer_a_dict = {}
    current_layer_a_dict = {}
    cur_depth = 0
    done_bound_list = []
    curr_best = -1
    # fixed values if Z is fixed
    # n = [] # total number of elements in each z strata
    # n1 = [] # total number of elements in each z strata for which outcome \in outcome_vals
    
    def __init__(self, data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, per_strata_a_counts_ch = [], safe_calculus = True, include_unmutated_treats=False):
        super().__init__(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, include_unmutated_treats = include_unmutated_treats)
        ### PAY ATTENTION: at the init phase and on bound() funct, self.per_strata_a_counts_ch contains the strata infos of the parent
        # and after the call objective() it contains the strata infos of the current node (to be passed to the children)
        self.per_strata_a_counts_ch = per_strata_a_counts_ch
        self.per_strata_a_counts_pa = per_strata_a_counts_ch
        if not hasattr(self,"df"):
            self.df = pd.read_csv(self.data_path)
        # n and n1 are constant per  z strata, they are not affected by treatments
        ns_per_strata = EffectCalculator.get_z_strata_tot_values(self.df, self.cfnds, [], [], self.outcome, self.outcome_vals)
        OurBnB.n  = [el[1] for el in ns_per_strata]
        OurBnB.n1 = [el[0] for el in ns_per_strata]
        self.current_depth = 0
        self.safe_calculus = safe_calculus

    def __repr__(self) -> str:
        s = super().__repr__()
        s += "self.per_strata_a_counts_pa:  "+str(self.per_strata_a_counts_pa)
        if hasattr(self,"per_strata_a_counts_ch"):
            s += "\nself.per_strata_a_counts_ch:  "+str(self.per_strata_a_counts_ch)
        if hasattr(self,"per_strata_a_counts_opt"):
            s += "\nself.per_strata_a_counts_opt: "+str(self.per_strata_a_counts_opt )
        return s

    #
    # required methods
    #
    def objective(self):        
        self.df = pd.read_csv(self.data_path)
        if len(self.treats) == 0:
            return -math.inf
        
        # must keep the Laplace correction to work with the tight optimistic estimator
        # after the call self.per_strata_a_counts contains the value of this node's strata and must be passed to child for optimistic bound estimation
        effect, self.per_strata_a_counts_ch = EffectCalculator.calculate_effect(self.df, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, laplace_correction = True, return_per_strata_a_counts = True)
        if len(self.treats) >= 2:
            if (np.array(self.per_strata_a_counts_opt) - np.array(self.per_strata_a_counts_ch) < 0).any():  #means I have more a in children that parents

                prev_lv_keys = [str([vv for vv in self.treats if vv!=v]) for v in self.treats]
                prev_lev_a = [[] for i in range(len(self.per_strata_a_counts_pa))]
                for k in prev_lv_keys:
                    dd = OurBnB.previous_layer_a_dict.get(k,[])
                    for i in range(len(self.per_strata_a_counts_pa)):
                        if dd == []:
                            #val = math.inf
                            val = 0   # if I don't have it calculated, since I'm running a breadth first, it means I discarded the rule
                        else:
                            val = dd[i]
                        prev_lev_a[i].append(val)

        OurBnB.current_layer_a_dict[str(self.treats)+str(self.treats_vals)] = self.per_strata_a_counts_ch  # adding a to the dict of saved As for tighter bound estimation
        OurBnB.previous_layer_a_dict[str(self.treats)+str(self.treats_vals)] = self.per_strata_a_counts_ch  # adding a to the dict of saved As for tighter bound estimation
        return effect
    
    def bound(self):
        OurBnB.done_bound_list.append(str(self.treats)+str(self.treats_vals))
        if len(self.treats) < 2:  # I don't have a value set up for per_strata_a_counts
            return self.unbounded_objective()
        optimist_effect = 0
        tot_cnt = 0

        prev_lv_keys = [str([self.treats[idx_i] for idx_i in range(len(self.treats)) if idx_i!=idx_j])+
            str([self.treats_vals[idx_i] for idx_i in range(len(self.treats)) if idx_i!=idx_j])
            for idx_j in range(len(self.treats))]
            
        prev_lev_a = [[] for i in range(len(self.per_strata_a_counts_pa))]
        for k in prev_lv_keys:
            dd = OurBnB.previous_layer_a_dict.get(k,[])
            for i in range(len(self.per_strata_a_counts_pa)):
                if dd == []:
                    if k in OurBnB.done_bound_list: # parent may be under calculus
                        val = math.inf
                    else:  # parent got pruned
                        if self.safe_calculus:
                            val = math.inf
                        else:
                            val = 0   # if I don't have it calculated, since I'm running a breadth first, it means I discarded the rule
                else:
                    val = dd[i]
                prev_lev_a[i].append(val)
        self.per_strata_a_counts_opt = [int(np.min([self.per_strata_a_counts_pa[i]]+prev_lev_a[i])) for i in range(len(self.per_strata_a_counts_pa))]

        # when I call bound() self.per_strata_a_counts_opt contains the value of each father's strata
        for i in range(len(self.per_strata_a_counts_opt)):
            a = self.per_strata_a_counts_opt[i]
            strata_best_eff = - math.inf
            for aa in range(a+1):
                curr_eff = (aa+1)/(aa + 2) - (OurBnB.n1[i]-aa+1)/(OurBnB.n[i]-aa+2) # this implicitely assumes Laplace correction
                if curr_eff > strata_best_eff:
                    strata_best_eff = curr_eff
            optimist_effect += strata_best_eff*OurBnB.n[i]
            tot_cnt += OurBnB.n[i]
        optimist_effect = optimist_effect/tot_cnt

        return optimist_effect
    
    def save_state(self, node):
        node.state = (self.data_path, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, self.per_strata_a_counts_pa, self.current_depth)

    def load_state(self, node):
        (self.data_path, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, self.per_strata_a_counts_pa, self.current_depth) = node.state

    def branch(self):
        if len(self.treats) + max(1,len(self.cfnds)) < OurBnB.max_rule_len: # with the max I consider the emptyset as a set of dim. 1 (since we correct on all single Zs and the empty set by construction)
            for v in OurBnB.total_treats:
                if v not in self.treats:
                    v_vals_list = [1]
                    if self.include_unmutated_treats:
                        v_vals_list.append(0)
                    for v_val in v_vals_list:
                        child = pybnb.Node()
                        if self.treats == []:
                            child_treats = [v]
                            child_treats_vals = [[v_val]]
                        else:
                            child_treats = self.treats + [v]
                            child_treats_vals = self.treats_vals + [[v_val]]
                        child.state = (self.data_path, self.cfnds, child_treats, child_treats_vals, self.outcome, self.outcome_vals, self.per_strata_a_counts_ch, self.current_depth + 1)
                        yield child


class ReliablePaperBnB(Simple):
    
    def __init__(self, data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, alpha, per_strata_data = [],
        use_graph_per_children = False, G = None, include_unmutated_treats = False):
        super().__init__(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len, include_unmutated_treats= include_unmutated_treats)
        ### PAY ATTENTION: at the init phase and on bound() funct, self.per_strata_data contains the strata infos of the parent
        # and after the call objective() it contains the strata infos of the current node (to be passed to the children)
        self.per_strata_data = per_strata_data
        self.alpha = alpha
        self.use_graph_per_children = use_graph_per_children
        self.G = G
            
    def __repr__(self) -> str:
        s = super().__repr__()
        s += "len(self.per_strata_data):"+str(len(self.per_strata_data))
        return s

    #
    # required methods
    #
    
    def objective(self):
        self.df = pd.read_csv(self.data_path)
        if len(self.treats) == 0:
            return -math.inf
        
        # must keep the Laplace correction to work with the tight optimistic estimator
        # after the call self.per_strata_data contains the value of this node's strata and must be passed to child for optimistic bound estimation
        effect, self.per_strata_data = EffectCalculator.calculate_effect(self.df, self.cfnds, 
            self.treats, self.treats_vals, self.outcome, self.outcome_vals, laplace_correction = True, 
            return_per_strata_data = True, apply_reliable_correction = True, correction_alpha = self.alpha)
        return effect
    
    def bound(self):
        # print(self.treats,self.treats_vals)

        if len(self.treats) < 2 or self.per_strata_data==[]:  # I don't have a value set up for per_strata_data
            return math.inf 
        optimist_effect = 0
        tot_cnt = 0
        # when I call bound() self.per_strata_data contains the value of each father's strata
        for pair in self.per_strata_data:
            a,b,c,d,n = pair
            n1 = a + c
            strata_best_eff = - math.inf
            for aa in range(a+1):
                curr_eff = (aa+1)/(aa + 2) - (n1-aa+1)/(n-aa+2) # this implicitely assumes Laplace correction
                curr_eff -= EffectCalculator.apply_correction(self.alpha, 2*math.sqrt(aa+2), 2*math.sqrt(n-aa+2))

                if curr_eff > strata_best_eff:
                    strata_best_eff = curr_eff
            optimist_effect += strata_best_eff*n
            tot_cnt += n
        optimist_effect = optimist_effect/tot_cnt

        return optimist_effect
    
    def save_state(self, node):
        node.state = (self.data_path, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, self.per_strata_data)

    def load_state(self, node):
        (self.data_path, self.cfnds, self.treats, self.treats_vals, self.outcome, self.outcome_vals, self.per_strata_data) = node.state

    def branch(self):
        v_vals_list = [1]
        if self.include_unmutated_treats:
            v_vals_list.append(0)

        if len(self.treats) + max(1,len(self.cfnds)) < ReliablePaperBnB.max_rule_len: # with the max I consider the emptyset as a set of dim. 1 (since we correct on all single Zs and the empty set by construction)
            if self.use_graph_per_children and len(self.treats)>0:

                # adj = [[n+"_somatic" for n in self.G.neighbors(t.split("_somatic")[0]) if n+"_somatic" not in self.treats 
                #         and n+"_somatic" not in self.cfnds] for t in self.treats if self.G.has_node(t.split("_somatic")[0])]
                # adj = [item for sublist in adj for item in sublist]

                adj = graph_utils.get_children_w_graph(treats=self.treats, total_treats=self.total_treats,
                    G = self.G, Zs=self.cfnds, data_cols=self.df.columns)
                for v in adj:
                    if v in self.df.columns and ReliablePaperBnB.total_treats.index(v)>max([ReliablePaperBnB.total_treats.index(vv) for vv in self.treats]+[-1]):
                        for v_val in v_vals_list:
                            child = pybnb.Node()
                            child_treats = self.treats + [v]
                            child_treats_vals = self.treats_vals + [[v_val]]
                            child.state = (self.data_path, self.cfnds, child_treats, child_treats_vals, self.outcome, self.outcome_vals, self.per_strata_data)
                            yield child

            else:
                for v in ReliablePaperBnB.total_treats:
                    if v not in self.treats and ReliablePaperBnB.total_treats.index(v)>max([ReliablePaperBnB.total_treats.index(vv) for vv in self.treats]+[-1]):
                        for v_val in v_vals_list:
                            child = pybnb.Node()
                            if self.treats == []:
                                child_treats = [v]
                                child_treats_vals = [[v_val]]
                            else:
                                child_treats = self.treats + [v]
                                child_treats_vals = self.treats_vals + [[v_val]]
                            child.state = (self.data_path, self.cfnds, child_treats, child_treats_vals, self.outcome, self.outcome_vals, self.per_strata_data)
                            yield child
    

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Select the class passing as additional argument 'OurBnB', 'OurBnB_single' or 'PaperBnB' and pass the max_rule_len value")
    data_path = "csvs/paper_example_rebranded_ext.csv"
    cfnds = ["Z"] 
    treats = []
    treats_vals = [[]]
    total_treats = ["E"+str(i+1) for i in range(40)] + ["X"]
    # total_treats = ["E1","E2","E3","E6","E7"] + ["X"]
    outcome = 'Y'
    outcome_vals = [5]
    # outcome_vals = [4]
    max_rule_len = int(sys.argv[2])

    if sys.argv[1] == 'OurBnB':
        try:
            problem = OurBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len)
            solver = pybnb.Solver()
            # solver = my_solver.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5, queue_strategy = "breadth"
                                , time_limit = 1000 # = 10 seconds
                                )
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)
        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))

    elif sys.argv[1] == 'OurBnB_single':
        try:
            problem = OurBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len, safe_calculus=False)
            solver = pybnb.Solver()
            # solver = my_solver.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5, queue_strategy = "breadth"
                                , time_limit = 1000 # = 10 seconds
                                )
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)

        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))
    
    elif sys.argv[1] == 'PaperBnB':
        try:
            problem = PaperBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len)
            solver = pybnb.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5)
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)

        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))
    
    elif sys.argv[1] == 'ReliablePaperBnB':
        try:
            problem = PaperBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len)
            solver = pybnb.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5)
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)

        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))

        try:
            problem = ReliablePaperBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len, alpha = 0.05, 
                include_unmutated_treats=True)
            solver = pybnb.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5)
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)

        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))
    
    elif sys.argv[1] == 'FigliaConIGraphi':
        data_path = "csvs/clinical_germline_somatictop515_molecular_subtype.csv"
        outcome = 'molecular_subtype'
        outcome_vals = [1]
        max_rule_len = int(sys.argv[2])

        split_column = 27
        data = pd.read_csv(data_path)

        # Select column values
        data_colnames = data.columns.values.tolist()
        # Isolate covariates
        Zs = "gender"
        # Drop last column as it is the outcome 
        data_colnames = data_colnames[:-1]
        # Isolate treatments
        Xs = data_colnames[split_column:]
        #y_vals = set(data[Y_variable])

        cfnds = [Zs]

        treats = []
        treats_vals = [[]]
        # treats = ["X3", "X5"]
        # treats_vals = [[1],[6]]
        # X3 \in [1] AND X5 \in [6] -> outcome \in outcome_vals | Z
        total_treats = Xs

        try:
            problem = ReliablePaperBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len,
                alpha=2.68791192e-9, use_graph_per_children=True)
            solver = pybnb.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5)
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)

        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))
    
    elif sys.argv[1] == 'TestGrafi0s':
        data_path =  "inputs/aggregated_mut_3N.csv"
        outcome = 'triple_negative'
        outcome_vals = [1]
        max_rule_len = int(sys.argv[2])

        split_column = 6
        data = pd.read_csv(data_path)

        # Select column values
        data_colnames = data.columns.values.tolist()
        # Isolate covariates
        Zs = "gender"
        # Drop last column as it is the outcome 
        data_colnames = data_colnames[:-1]
        # Isolate treatments
        Xs = data_colnames[split_column:]
        #y_vals = set(data[Y_variable])

        cfnds = [Zs]

        treats = []
        treats_vals = [[]]
        # treats = ["X3", "X5"]
        # treats_vals = [[1],[6]]
        # X3 \in [1] AND X5 \in [6] -> outcome \in outcome_vals | Z
        total_treats = Xs

        df_graph = pd.read_table('PPI/FIsInGene_122921_with_annotations.txt')
        G = nx.from_pandas_edgelist(df_graph, 'Gene1', 'Gene2', 'Score')

        try:
            problem = ReliablePaperBnB(data_path, cfnds, treats, treats_vals, outcome, outcome_vals, total_treats, max_rule_len = max_rule_len,
                alpha=2.68791192e-9, use_graph_per_children=True, G = G, include_unmutated_treats = True)
                                        
            solver = pybnb.Solver()
            results = solver.solve(problem,
                                absolute_gap=1e-5)
            # print(results)
            if solver.is_dispatcher:
                print(results.best_node.state)

        except OverflowError:
            print("Overflow Error TOT", sys.getsizeof(Simple.SAVED_STATE))
    
    else:
        print("Wrong type")
    