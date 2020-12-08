import numpy as np
import xml.etree.ElementTree as ET
try:
	__IPYTHON__
except NameError:
	in_notebook = False
else:
	in_notebook = True

""" test for tqdm progress bars """

try:
	if in_notebook:
		from tqdm import tqdm_notebook as tqdm
	else:
		from tqdm import tqdm
	has_tqdm=True
except:
	has_tqdm=False



GlobalGroup = set(range(48)) # Assume Oh

class Vertex:
	def __init__(self):
		self.E = 0.0
		self.position = np.zeros(3)
		self.ku = 0.0
		self.time = 1.0
		self.temperature = 0.0
		self.TADBarrier = 0.0
		self.label = np.uint64(0)
		self.clabel = np.uint64(0)
		self.selfsyms=[]
		self.conjugate=[]
		self.equivalents=[]
	def shift(self,pos):
		self.position -= pos.copy()

	def status(self,pbc=None,verbose=False):
		res = "\nLabel: %s %s\n" % (str(self.clabel),str(self.label))
		res += "\tE: %f\n" % self.E
		res += "\tku: %10.10g\n\n" % (self.ku)
		correc = np.zeros(3)

		if not pbc is None:
			po = pbc(self.position,c=0.0)
		else:
			po = self.position
		res += "pos: %s\n\n" % np.array2string(po)

		for ops in self.selfsyms:
			M=opMat(ops[0])
			#for i in range(3):
			#    res += "\t"+np.array2string(M[i])+"\n"
			if not pbc is None:
				ss = pbc(np.r_[ops[1:]])
				dd = pbc( opMat(ops[0]).dot(self.position)-self.position+ss ,c=1.0)
			else:
				ss = np.r_[ops[1:]]
				dd = (opMat(ops[0])).dot(self.position)-self.position+ss
			res += "M:\n%s,\tpbc[shi]:%s\t,\tpbc[(M-I).pos+shi] = %s\n" % (np.array2string(M),np.array2string(ss),np.array2string(dd))


		#for ops in self.equivalents:
		#    res += "\tSeenEquivalents: %s, %d, [%f %f %f]\n" % (ops[0],ops[1],ops[2],ops[3],ops[4])
		res += "\n\n\n"
		return res

class Edge:
	def __init__(self):
		self.iclabel=np.uint64(0)
		self.fclabel=np.uint64(0)
		self.ilabel=np.uint64(0)
		self.flabel=np.uint64(0)
		self.fbar=0.0
		self.bbar=0.0
		self.fnu=0.0
		self.bnu=0.0
		self.selftran=None
		self.equivalents=[]

	def status(self,center=None):
		res = "\nLabels: (%s,%s) -> (%s,%s)\n" % (str(self.iclabel),str(self.ilabel),str(self.fclabel),str(self.flabel))
		res += "\tSelf Transition? %d\n" % bool(not self.selftran is None)
		res += "\tdE: %f %f\n" % (self.fbar,self.bbar)
		if not self.selftran is None:
			op = self.selftran[0]
			v = np.r_[[self.selftran[1],self.selftran[2],self.selftran[3]]]
			if not center is None:
				v -= center
			res += "\tJump Transform: %d, [%f %f %f]\n" % (op,v[0],v[1],v[2])

		if len(self.equivalents) > 0:
			res += "\tFound Equivalents: \n"
			for ops in self.equivalents:
				res += "\t%s -> %s, %d, [%f %f %f]\n" % (ops[0],ops[1],ops[2],ops[3],ops[4],ops[5])
		res += "\n"
		return res

class DiffusionModel:
	def __init__(self,state_file='data/MarkovModel.xml',correct=False):
		"""
			Read in XML file which specifies:
				- Primitive Unit Cell (PUC)
				- Conventional Unit Cell (CUC)
				- Simulation Supercell in units or matrix multiplication of CUC
		"""
		self.PrimitiveUnitCell = None
		self.UnitCell = None
		self.SuperCell = None
		self.PointGroup = GlobalGroup

		""" parse xml file """
		self.Cell = np.zeros((3,3))
		self.center = np.zeros(3)
		self.states = []
		self.transitions = []

		self.LatticeConstant = -1.0

		self.nstates=0
		self.index = {}
		self.cmap = {}

		tree = ET.parse(state_file)
		for child in tree.getroot():
			if child.tag=='Lattice':
				data = child.text.strip().split(" ")
				self.LatticeConstant = float(data[0])
				print("\n\tLattice Constant = %fA" % self.LatticeConstant)
				if len(data)>1:
					if data[1] == 'bcc':
						self.PrimitiveUnitCell = 0.5*np.ones((3,3))-np.eye(3)
					if data[1] == 'fcc':
						self.PrimitiveUnitCell = 0.5*np.ones((3,3)) - 0.5*np.eye(3)
					print("\n\tLattice Structure = %s" % data[1])

			if child.tag=='PrimitiveUnitCell' and (self.PrimitiveUnitCell is None):
				self.PrimitiveUnitCell = np.genfromtxt(child.text.splitlines(),dtype=np.float)
			if child.tag=='UnitCell' and (self.UnitCell is None):
				self.UnitCell = np.genfromtxt(child.text.splitlines(),dtype=np.float)

			if child.tag=='SuperCell':
				self.SuperCell = np.genfromtxt(child.text.splitlines(),dtype=np.float)

			if not ((self.PrimitiveUnitCell is None) or (self.UnitCell is None) or (self.SuperCell is None) or (self.LatticeConstant<0.0)):
				self.PrimitiveUnitCell *= self.LatticeConstant
				self.UnitCell *= self.LatticeConstant
				break

		if self.SuperCell.size==3:
			self.Cell = self.UnitCell.dot(np.diag(self.SuperCell))
		else:
			self.Cell = self.SuperCell.dot(self.PrimitiveUnitCell)

		self.center = self.Cell.dot(np.ones(3))/2.0

		for child in tree.getroot():
			if child.tag=='Vertex':
				ss = Vertex()
				for ele in child:
					if ele.tag == 'CanonLabel':
						ss.clabel = np.uint64(ele.text)
					elif ele.tag == 'ReferenceLabel':
						ss.label = np.uint64(ele.text)
					elif ele.tag == 'Energy':
						ss.E = np.float64(ele.text)
					elif ele.tag == 'Position':
						dd = np.fromstring(ele.text,dtype=np.float64,sep=' ')
						ss.position = dd.copy()
					elif ele.tag == 'Time':
						ss.time = np.float64(ele.text)
					elif ele.tag == 'UnknownRate':
						ss.ku = np.float64(ele.text)
					elif ele.tag == 'Temperature':
						ss.temperature = np.float64(ele.text)
					elif ele.tag == 'TADBarrier':
						ss.TADBarrier = np.float64(ele.text)
					elif ele.tag == 'SelfSymmetries':
						dd = np.genfromtxt(ele.text.splitlines(),dtype=(np.int,np.float,np.float,np.float))
						#np.fromstring(ele.text,dtype=np.float64,sep=' ')
						#dd[1:] = self.pucsnap(dd[1:])
						for _dd in dd:
							ss.selfsyms.append([_dd[0],_dd[1],_dd[2],_dd[3]])
					elif ele.tag == 'SeenEquivalents':
						sid = np.genfromtxt(ele.text.splitlines(),dtype=np.uint64)
						dd = np.genfromtxt(ele.text.splitlines())

						if dd.size==dd.shape[0]:
							dd = dd.reshape((1,-1))
							sid = sid.reshape((1,-1))

						#,dtype=(np.uint64,np.int,np.float,np.float,np.float))#,np.float,np.float,np.float))
						for _dd in zip(sid,dd):
							ss.equivalents.append([_dd[0][0].astype(np.uint64),_dd[0][1].astype(np.int),_dd[1][2],_dd[1][3],_dd[1][4]])
							#list(_dd))
				ss.conjugate =  self.conjugate_set(ss).copy()
				self.states.append(ss)
			elif child.tag=='Edge':
				ee = Edge()
				for ele in child:
					if ele.tag == 'CanonLabels':
						dd = np.fromstring(ele.text,dtype=np.uint64,sep=' ')
						ee.iclabel,ee.fclabel = dd[0],dd[1]
					if ele.tag == 'ReferenceLabels':
						dd = np.fromstring(ele.text,dtype=np.uint64,sep=' ')
						ee.ilabel,ee.flabel = dd[0],dd[1]
					elif ele.tag == 'Barriers':
						dd = np.fromstring(ele.text,dtype=np.float64,sep=' ')
						ee.fbar, ee.bbar = dd[0],dd[1]
					elif ele.tag == 'PreFactors':
						dd = np.fromstring(ele.text,dtype=np.float64,sep=' ')
						ee.fnu, ee.bnu = dd[0],dd[1]
					elif ele.tag == 'TransitionSymmetry':
						dd = np.fromstring(ele.text,dtype=np.float64,sep=' ')
						ee.selftran = [int(dd[0]),dd[1],dd[2],dd[3]]
					elif ele.tag == 'SeenEquivalents':
						lli = np.uint64(ele.text.strip().split(" ")[0])
						llf = np.uint64(ele.text.strip().split(" ")[1])
						dd = np.fromstring(ele.text,dtype=np.float64,sep=' ')[2:]
						ee.equivalents.append([lli,llf,int(dd[0]),dd[1],dd[2],dd[3]])
				self.transitions.append(ee)
		print("\n\tPrimitive Unit Cell:",end="")
		for ii in range(3):
			print("\n\t",self.PrimitiveUnitCell[ii])
		print("\n\tSuperCell:\n\t",self.Cell[0],"\n\t",self.Cell[1],"\n\t",self.Cell[2])
		print("\n\tCenter:\n\t",self.center,"\n")

		if correct:
			iC = np.linalg.inv(self.PrimitiveUnitCell)
			C = self.PrimitiveUnitCell
			ms=0.0
			com = np.zeros(3)
			for i in range(len(self.states)):
				nsym = len(self.states[i].selfsyms)
				dp = np.zeros(3)
				for j in range(nsym):
					M = opMat(self.states[i].selfsyms[j][0])
					c = C.dot(np.round(2.0*iC.dot(self.states[i].selfsyms[j][1:]))/2.0)
					op = self.pbcuc(self.states[i].position.copy(),c=0.0)
					tp = self.pbcuc(M.dot(op)+c,c=0.0)
					dp += (tp-op)/float(nsym)
				self.states[i].position += dp
		self.fill_cmap_index()

	def conjugate_set(self,s):
		r"""
			Find nonunique conjugate set :math:`C` to coset :math:`S` for each element
			such that :math:`{C}\times{S}` gives the lattice point group :math:`G`

			The conjugate set thus enumerates all nonequivalent instances of a state

			We arbitrarily pick left cosets (must be consistent throughout)

			I.e. a state :math:`A` is invariant under any element :\math:`s\in{S}`
			such that :math:`sA\equiv{A}`.

			For any element :math:`g\in{G}` we can always make a decomposition
			:math:`g={c}s` but we must find :math:`c` for each :math:`g\in{G}`.

			So we want set :math:`C` that gives :math:`G` when *pre*multiplying
			every member of :math:`S` by every member of :math:`C`, which is what
			this function finds
		"""
		S = set([ss[0] for ss in s.selfsyms]) ## selfsymmetry of state
		P, N = set([0]), set([0])
		# While Union(P,S) does not have all the elements of G
		while not (P|S)>=self.PointGroup:
			# First element "op" in G but not P or S
			for op in (self.PointGroup-P)-S:
				# ss = op x S
				ss = set([opMult(op,s) for s in S])
				# Add to P all elements in ss that are not in S
				P |= ss-S # Ensures intersection(P,S) = [0]
				# Add op to N
				N |= set([op])
				break # one at a time....
		return N

	def transform_index(self,state,operation):
		r"""
			When acting on a reference state :math:`A` with :math:`g\in{G}`,
			which nonequivalent state does it return?

			If :math:`g\in{S}`, i.e. in the (left) coset of :\math:`A`, then
			we recover the same state. In general, however, this will not be the case

			So we have :math:`gA=csA=cA\forall s\in{S}`

			The intersection of the set :math:`c\prime=\{gs^{-1} \forall s\in{S}\}`
			with the conjugate set :math:`C` should therefore return one element
		"""

		ss = set([j[0] for j in self.states[state].selfsyms]) # S
		cj = set(self.states[state].conjugate.copy()) # C
		nn = set(opMult(operation,opInv(s)) for s in ss) # op s^{-1}
		#nn = set(opMult(opInv(s),operation) for s in ss)

		if len(nn&cj)!=1:
			print("DEGEN!!",len(nn&cj))
			return -1
		return list(nn&cj)[0]

	def fill_cmap_index(self):
		r"""
			Establish index for each nonequivalent state from group notation
		"""
		self.nstates=0
		self.index = {}
		self.cmap = {}
		iC = np.linalg.inv(self.UnitCell)
		C = self.UnitCell
		for i,s in enumerate(self.states):
			imap={}
			imap[s.label] = [0,s.position] # np.zeros(3) ]
			for ee in s.equivalents:
				op = self.transform_index(i,ee[1])
				dp = self.pbc(np.r_[ee[-3:]]-s.position,c=1.0)
				dp = 0.5*C@np.round(2.0*iC@dp)
				imap[ee[0]] = [op,dp,ee[1]]
			self.cmap[s.clabel] = [i,imap.copy()]
		for i,s in enumerate(self.states):
			lind = {}
			for c in s.conjugate:
				lind[c] = self.nstates
				self.nstates += 1
			self.index[i] = lind
		print("\n\tDecompressed system has %d states irredudible under translation" % self.nstates)


	def fill_rates(self,evalT):
		"""
			Evaluate rate matrix K and displacement tensors Kxd, Kxdxd
		"""
		beta = 1.0/8.617e-5/evalT
		iC = np.linalg.inv(self.UnitCell)
		C = self.UnitCell
		pi = np.zeros(self.nstates)
		Ku = np.zeros(self.nstates)
		K = np.zeros((self.nstates,self.nstates))
		Kd = np.zeros((self.nstates,self.nstates,3))
		Kdd = np.zeros((self.nstates,self.nstates,3,3))

		for ss in self.states:
			icl= self.cmap[ss.clabel][0]
			tad_boltz = ss.TADBarrier / 8.617e-5 / ss.temperature
			ku = ss.ku * np.exp(-(ss.temperature/evalT-1.0)*tad_boltz)
			for j in ss.conjugate:
				i = self.index[icl][j]
				Ku[i] = ku
				pi[i] = np.exp(-beta*ss.E)
		pi/=pi.sum()


		for ee in self.transitions:

			icl = self.cmap[ee.iclabel][0]
			il  = self.cmap[ee.iclabel][1][ee.ilabel][0]
			p   = self.cmap[ee.iclabel][1][ee.ilabel][1]
			fcl = self.cmap[ee.fclabel][0]
			fl  = self.cmap[ee.fclabel][1][ee.flabel][0]
			fp  = self.cmap[ee.fclabel][1][ee.flabel][1]
			dp  = self.pbc(fp-p,c=1.0)
			rdp = C.dot(np.round(2.0*iC.dot(dp)))/2.0
			kf = ee.fnu*np.exp(-ee.fbar*beta)
			kb = ee.bnu*np.exp(-ee.bbar*beta)

			for sop in self.PointGroup:
				i = self.index[icl][self.transform_index(icl,opMult(sop,il))]
				f = self.index[fcl][self.transform_index(fcl,opMult(sop,fl))]
				d = opMat(sop).dot(rdp)

				K[f][i] += kf
				Kd[f][i] += d*kf
				Kdd[f][i] += kf*np.outer(d,d)

				K[i][f] += kb
				Kd[i][f] -= d*kb
				Kdd[i][f] += kb*np.outer(d,d)

		return K,Kd,Kdd,Ku,pi

	def calcD(self,nK,nKd,nKdd,nku,pi,ls=True,rn=False,startP=None):
		"""
			Worker function to evaluate diffusion matrix from rate tensors;
			front end functions are calculateD and sampleD
		"""

		iG = np.diag(nK.sum(axis=0)+nku)-nK

		if ls:
			x = np.linalg.solve(iG,pi)
		else:
			x = np.linalg.pinv(iG).dot(pi)

		tau = x.sum()

		if rn:
			dallas_mu = np.linalg.solve(iG.T,nKd.sum(axis=0))/tau
			nv = np.diag(nK.sum(axis=0)+nku)@x
			nv /= nv.sum()
		x /= x.sum()


		if not startP is None:
			x = startP
		Dc = np.zeros((3,3))
		mu = np.tensordot(nKd.sum(axis=0),x,axes=(0,0))
		GKdP = np.linalg.solve(iG,np.tensordot(nKd,x,axes=(1,0)))
		Dc = np.tensordot(nKd.sum(axis=0),GKdP,axes=(0,0)) - np.outer(mu,mu)
		Duc = np.tensordot(nKdd.sum(axis=0),x,axes=(0,0))
		if rn:
			return Duc, Dc+Dc.T, tau, dallas_mu, mu,nv
		else:
			return Duc, Dc+Dc.T, tau


	# Ks = K.PI
	def sampleK(self,eK,eKd,eKdd,eKu,epi):
		eN = eKu.size
		nK = eK.copy()
		nKd = eKd.copy()
		nKdd = eKdd.copy()
		npi = epi.copy()
		kup = eKu*epi
		ll = np.tril_indices(eN)
		dKs = np.zeros((eN,eN))

		for i in np.argsort(np.random.uniform(size=len(ll[0]))):
			l = [ll[0][i],ll[1][i]]
			dkup = (kup-dKs.sum(axis=0))*0.95
			dKs[l[0],l[1]] = np.random.uniform(0.95)*min(dkup[l[0]],dkup[l[1]])
			dKs[l[1],l[0]] = dKs[l[0],l[1]]

		dK = dKs.dot(np.diag(1.0/npi))
		nKu = eKu.copy()-dK.sum(axis=0)
		nK += dK.copy()

		for l in zip(ll[0],ll[1]):
			di = np.random.randint(4)
			if di<3:
				sdv = self.PrimitiveUnitCell[:,di]
			else:
				sdv = np.zeros(3)
			sdvv = np.outer(sdv,sdv)
			nKd[l[0],l[1]] += dK[l[0],l[1]]*sdv
			nKd[l[1],l[0]] -= dK[l[1],l[0]]*sdv
			nKdd[l[0],l[1]] += dK[l[0],l[1]]*sdvv
			nKdd[l[1],l[0]] += dK[l[1],l[0]]*sdvv
		Duc,Dc,tau = self.calcD(nK,nKd,nKdd,nKu,npi,ls=True)
		return Duc+Dc,tau

	def calculateD(self,evalT=100.0,ev=True,pure=False):
		"""
			Evaluate diffusion matrix and residence time at a given temperature
		"""
		K,Kd,Kdd,Ku,pi = self.fill_rates(evalT)

		if pure: # simulate ku->0 limit
			Ku *= 1.0e-10
		Duc,Dc,tau = self.calcD(K,Kd,Kdd,Ku,pi,ls=True)
		if ev:
			return np.linalg.eigvalsh(Duc+Dc),tau
		else:
			return Duc+Dc,tau


	def sampleD(self,evalT=100.0,Nsample=10):
		"""
			Returns Nsample projections of modified diffusion tensors
			onto the eigenvectors of the original diffusion tensor
			When the modifed and original diffusion tensor can be
			simultaneously diagonalized this amounts an extraction of modified
			eigenvalues
		"""
		DIM=3
		beta = 1.0/evalT/8.617e-5

		tau_mmm = np.zeros(3) # mean max min
		D_mmm = np.zeros((3,3)) # mean max min\s

		Dsamarr = np.zeros((Nsample,DIM+1))

		K,Kd,Kdd,Ku,pi = self.fill_rates(evalT)
		Duc,Dc,tau_mmm[0] = self.calcD(K,Kd,Kdd,Ku,pi,ls=True)
		D_mmm[0],Dvec = np.linalg.eigh(Duc+Dc)
		if has_tqdm:
			pbar = tqdm(total=Nsample,leave=False,desc="MC-UQ, T=%dK" % evalT)
		for i in range(Nsample):
			Ds,Dsamarr[i][-1] = self.sampleK(K,Kd,Kdd,Ku,pi)
			Dsamarr[i][:-1] = np.diag((Dvec.T)@Ds@Dvec)
			if has_tqdm:
				pbar.update(1)
		if has_tqdm:
			pbar.close()

		# Take largest sum for tau
		tau_mmm[1] = Dsamarr[:,3].max()
		tau_mmm[2] = Dsamarr[:,3].min()
		D_mmm[1] = Dsamarr[:,:3].max(axis=0)
		D_mmm[2] = Dsamarr[:,:3].min(axis=0)

		return D_mmm,tau_mmm,Dvec



	def pbc(self,v,sel=None,c=True):
		if sel is None:
			return v-np.dot(self.Cell,np.floor(np.dot(np.linalg.inv(self.Cell),v.T)+0.5*c)).T
		else:
			return v-np.dot(self.Cell[:,sel][sel,:],np.floor(np.dot(np.linalg.inv(self.Cell[:,sel][sel,:]),v.T)+0.5*c)).T
	def pbcuc(self,v,c=False):
		return v-np.dot(self.PrimitiveUnitCell,np.floor(np.dot(np.linalg.inv(self.PrimitiveUnitCell),v.T)+0.5*c)).T

	def pbccuc(self,v,c=False):
		return v-np.dot(self.UnitCell,np.floor(np.dot(np.linalg.inv(self.UnitCell),v.T)+0.5*c)).T
	def pucsnap(self,v):
		return np.dot(self.PrimitiveUnitCell,0.5*np.round(2.0*np.dot(np.linalg.inv(self.PrimitiveUnitCell),v.T))).T

multtable = np.r_[[
	[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47],
	[1,2,0,5,3,4,7,8,6,11,9,10,13,14,12,17,15,16,19,20,18,23,21,22,25,26,24,29,27,28,31,32,30,35,33,34,37,38,36,41,39,40,43,44,42,47,45,46],
	[2,0,1,4,5,3,8,6,7,10,11,9,14,12,13,16,17,15,20,18,19,22,23,21,26,24,25,28,29,27,32,30,31,34,35,33,38,36,37,40,41,39,44,42,43,46,47,45],
	[3,4,5,0,1,2,9,10,11,6,7,8,15,16,17,12,13,14,21,22,23,18,19,20,27,28,29,24,25,26,33,34,35,30,31,32,39,40,41,36,37,38,45,46,47,42,43,44],
	[4,5,3,2,0,1,10,11,9,8,6,7,16,17,15,14,12,13,22,23,21,20,18,19,28,29,27,26,24,25,34,35,33,32,30,31,40,41,39,38,36,37,46,47,45,44,42,43],
	[5,3,4,1,2,0,11,9,10,7,8,6,17,15,16,13,14,12,23,21,22,19,20,18,29,27,28,25,26,24,35,33,34,31,32,30,41,39,40,37,38,36,47,45,46,43,44,42],
	[6,13,20,15,22,11,0,25,38,27,40,5,24,1,32,3,34,29,36,31,2,33,4,41,12,7,44,9,46,17,42,19,14,21,16,47,18,43,8,45,10,23,30,37,26,39,28,35],
	[7,14,18,17,21,10,1,26,36,29,39,4,25,2,30,5,33,28,37,32,0,35,3,40,13,8,42,11,45,16,43,20,12,23,15,46,19,44,6,47,9,22,31,38,24,41,27,34],
	[8,12,19,16,23,9,2,24,37,28,41,3,26,0,31,4,35,27,38,30,1,34,5,39,14,6,43,10,47,15,44,18,13,22,17,45,20,42,7,46,11,21,32,36,25,40,29,33],
	[9,16,23,12,19,8,3,28,41,24,37,2,27,4,35,0,31,26,39,34,5,30,1,38,15,10,47,6,43,14,45,22,17,18,13,44,21,46,11,42,7,20,33,40,29,36,25,32],
	[10,17,21,14,18,7,4,29,39,26,36,1,28,5,33,2,30,25,40,35,3,32,0,37,16,11,45,8,42,13,46,23,15,20,12,43,22,47,9,44,6,19,34,41,27,38,24,31],
	[11,15,22,13,20,6,5,27,40,25,38,0,29,3,34,1,32,24,41,33,4,31,2,36,17,9,46,7,44,12,47,21,16,19,14,42,23,45,10,43,8,18,35,39,28,37,26,30],
	[12,19,8,9,16,23,24,37,2,3,28,41,0,31,26,27,4,35,30,1,38,39,34,5,6,43,14,15,10,47,18,13,44,45,22,17,42,7,20,21,46,11,36,25,32,33,40,29],
	[13,20,6,11,15,22,25,38,0,5,27,40,1,32,24,29,3,34,31,2,36,41,33,4,7,44,12,17,9,46,19,14,42,47,21,16,43,8,18,23,45,10,37,26,30,35,39,28],
	[14,18,7,10,17,21,26,36,1,4,29,39,2,30,25,28,5,33,32,0,37,40,35,3,8,42,13,16,11,45,20,12,43,46,23,15,44,6,19,22,47,9,38,24,31,34,41,27],
	[15,22,11,6,13,20,27,40,5,0,25,38,3,34,29,24,1,32,33,4,41,36,31,2,9,46,17,12,7,44,21,16,47,42,19,14,45,10,23,18,43,8,39,28,35,30,37,26],
	[16,23,9,8,12,19,28,41,3,2,24,37,4,35,27,26,0,31,34,5,39,38,30,1,10,47,15,14,6,43,22,17,45,44,18,13,46,11,21,20,42,7,40,29,33,32,36,25],
	[17,21,10,7,14,18,29,39,4,1,26,36,5,33,28,25,2,30,35,3,40,37,32,0,11,45,16,13,8,42,23,15,46,43,20,12,47,9,22,19,44,6,41,27,34,31,38,24],
	[18,7,14,21,10,17,36,1,26,39,4,29,30,25,2,33,28,5,0,37,32,3,40,35,42,13,8,45,16,11,12,43,20,15,46,23,6,19,44,9,22,47,24,31,38,27,34,41],
	[19,8,12,23,9,16,37,2,24,41,3,28,31,26,0,35,27,4,1,38,30,5,39,34,43,14,6,47,15,10,13,44,18,17,45,22,7,20,42,11,21,46,25,32,36,29,33,40],
	[20,6,13,22,11,15,38,0,25,40,5,27,32,24,1,34,29,3,2,36,31,4,41,33,44,12,7,46,17,9,14,42,19,16,47,21,8,18,43,10,23,45,26,30,37,28,35,39],
	[21,10,17,18,7,14,39,4,29,36,1,26,33,28,5,30,25,2,3,40,35,0,37,32,45,16,11,42,13,8,15,46,23,12,43,20,9,22,47,6,19,44,27,34,41,24,31,38],
	[22,11,15,20,6,13,40,5,27,38,0,25,34,29,3,32,24,1,4,41,33,2,36,31,46,17,9,44,12,7,16,47,21,14,42,19,10,23,45,8,18,43,28,35,39,26,30,37],
	[23,9,16,19,8,12,41,3,28,37,2,24,35,27,4,31,26,0,5,39,34,1,38,30,47,15,10,43,14,6,17,45,22,13,44,18,11,21,46,7,20,42,29,33,40,25,32,36],
	[24,31,38,27,34,41,12,43,20,15,46,23,6,19,44,9,22,47,42,13,8,45,16,11,0,37,32,3,40,35,36,1,26,39,4,29,30,25,2,33,28,5,18,7,14,21,10,17],
	[25,32,36,29,33,40,13,44,18,17,45,22,7,20,42,11,21,46,43,14,6,47,15,10,1,38,30,5,39,34,37,2,24,41,3,28,31,26,0,35,27,4,19,8,12,23,9,16],
	[26,30,37,28,35,39,14,42,19,16,47,21,8,18,43,10,23,45,44,12,7,46,17,9,2,36,31,4,41,33,38,0,25,40,5,27,32,24,1,34,29,3,20,6,13,22,11,15],
	[27,34,41,24,31,38,15,46,23,12,43,20,9,22,47,6,19,44,45,16,11,42,13,8,3,40,35,0,37,32,39,4,29,36,1,26,33,28,5,30,25,2,21,10,17,18,7,14],
	[28,35,39,26,30,37,16,47,21,14,42,19,10,23,45,8,18,43,46,17,9,44,12,7,4,41,33,2,36,31,40,5,27,38,0,25,34,29,3,32,24,1,22,11,15,20,6,13],
	[29,33,40,25,32,36,17,45,22,13,44,18,11,21,46,7,20,42,47,15,10,43,14,6,5,39,34,1,38,30,41,3,28,37,2,24,35,27,4,31,26,0,23,9,16,19,8,12],
	[30,37,26,39,28,35,42,19,14,21,16,47,18,43,8,45,10,23,12,7,44,9,46,17,36,31,2,33,4,41,0,25,38,27,40,5,24,1,32,3,34,29,6,13,20,15,22,11],
	[31,38,24,41,27,34,43,20,12,23,15,46,19,44,6,47,9,22,13,8,42,11,45,16,37,32,0,35,3,40,1,26,36,29,39,4,25,2,30,5,33,28,7,14,18,17,21,10],
	[32,36,25,40,29,33,44,18,13,22,17,45,20,42,7,46,11,21,14,6,43,10,47,15,38,30,1,34,5,39,2,24,37,28,41,3,26,0,31,4,35,27,8,12,19,16,23,9],
	[33,40,29,36,25,32,45,22,17,18,13,44,21,46,11,42,7,20,15,10,47,6,43,14,39,34,5,30,1,38,3,28,41,24,37,2,27,4,35,0,31,26,9,16,23,12,19,8],
	[34,41,27,38,24,31,46,23,15,20,12,43,22,47,9,44,6,19,16,11,45,8,42,13,40,35,3,32,0,37,4,29,39,26,36,1,28,5,33,2,30,25,10,17,21,14,18,7],
	[35,39,28,37,26,30,47,21,16,19,14,42,23,45,10,43,8,18,17,9,46,7,44,12,41,33,4,31,2,36,5,27,40,25,38,0,29,3,34,1,32,24,11,15,22,13,20,6],
	[36,25,32,33,40,29,18,13,44,45,22,17,42,7,20,21,46,11,6,43,14,15,10,47,30,1,38,39,34,5,24,37,2,3,28,41,0,31,26,27,4,35,12,19,8,9,16,23],
	[37,26,30,35,39,28,19,14,42,47,21,16,43,8,18,23,45,10,7,44,12,17,9,46,31,2,36,41,33,4,25,38,0,5,27,40,1,32,24,29,3,34,13,20,6,11,15,22],
	[38,24,31,34,41,27,20,12,43,46,23,15,44,6,19,22,47,9,8,42,13,16,11,45,32,0,37,40,35,3,26,36,1,4,29,39,2,30,25,28,5,33,14,18,7,10,17,21],
	[39,28,35,30,37,26,21,16,47,42,19,14,45,10,23,18,43,8,9,46,17,12,7,44,33,4,41,36,31,2,27,40,5,0,25,38,3,34,29,24,1,32,15,22,11,6,13,20],
	[40,29,33,32,36,25,22,17,45,44,18,13,46,11,21,20,42,7,10,47,15,14,6,43,34,5,39,38,30,1,28,41,3,2,24,37,4,35,27,26,0,31,16,23,9,8,12,19],
	[41,27,34,31,38,24,23,15,46,43,20,12,47,9,22,19,44,6,11,45,16,13,8,42,35,3,40,37,32,0,29,39,4,1,26,36,5,33,28,25,2,30,17,21,10,7,14,18],
	[42,43,44,45,46,47,30,31,32,33,34,35,36,37,38,39,40,41,24,25,26,27,28,29,18,19,20,21,22,23,6,7,8,9,10,11,12,13,14,15,16,17,0,1,2,3,4,5],
	[43,44,42,47,45,46,31,32,30,35,33,34,37,38,36,41,39,40,25,26,24,29,27,28,19,20,18,23,21,22,7,8,6,11,9,10,13,14,12,17,15,16,1,2,0,5,3,4],
	[44,42,43,46,47,45,32,30,31,34,35,33,38,36,37,40,41,39,26,24,25,28,29,27,20,18,19,22,23,21,8,6,7,10,11,9,14,12,13,16,17,15,2,0,1,4,5,3],
	[45,46,47,42,43,44,33,34,35,30,31,32,39,40,41,36,37,38,27,28,29,24,25,26,21,22,23,18,19,20,9,10,11,6,7,8,15,16,17,12,13,14,3,4,5,0,1,2],
	[46,47,45,44,42,43,34,35,33,32,30,31,40,41,39,38,36,37,28,29,27,26,24,25,22,23,21,20,18,19,10,11,9,8,6,7,16,17,15,14,12,13,4,5,3,2,0,1],
	[47,45,46,43,44,42,35,33,34,31,32,30,41,39,40,37,38,36,29,27,28,25,26,24,23,21,22,19,20,18,11,9,10,7,8,6,17,15,16,13,14,12,5,3,4,1,2,0]]
	].astype(int)
"""
	std::array<double,NDIM*NDIM> M;
	for(int i=0;i<NDIM*NDIM;i++) fM[i] = 0.0;
	for(int i=0;i<NDIM*NDIM;i++) M[i] = 0.0;
	for(int i=0;i<NDIM;i++) M[i*NDIM+i] = 1.0;

	#if NDIM==3
	/*
		48 SYMMETRIES OF THE CUBE = 3!=6 permutations * 2^3=8 reflections
		operation = 6R + P,  R E [0,7],  P E [0,5]
		R = 0 nothing
		R = 1-3 single reflection x,y,z
		R = 4-6 double reflection xy yz zx NOTE ORDER
		R = 7 triple reflection
		P = 0 nothing
		P = 1-2 cyclic permutation: nothing, xyz -> yzx,  xyz -> zxy
		P = 3-5 x<->y then cyclic permutation
	*/
	if(op>=48) op=0;
	unsigned int R = op/6;
	unsigned int P = op%6;
	// Reflections
	// Diagonal elements 0,4,8
	// R = [1-3] : M[(R-1)*4] = -1 														==>    M[4*((R-1+i)%3)] = -1 for i = 0
	// R = [4-6] : M[(R-4)*4] = -1 and M[4*(R-4+1)%3] = -1    ==>    M[4*((R-1+i)%3)] = -1 for i = 0,1
	// R = 7 : M[0] = -1 and M[4] = -1 and M[8] = -1					==>    M[4*((R-1+i)%3)] = -1 for i = 0,1,2
	int ii = (R>0)*(1+(R-1)/3); // R->ii : 0->0, [1-3]->1, [4-6]->2, 7->3
	for(int i=0;i<ii; i++) M[4*((R-1+i)%3)] = -1.;
	// Permutations
	for(int i=0;i<3;i++) {
		// i->ii     P<3: 012->012     else:   012->102   (swap x<->y if P >= 3)
		// P<3:   P/3 =0, ii = i
		// P>=3:  P/3 =1, ii = 1-i%2 + i/2  [  i<2:   ii = 1-i  (01->10)      i==2:   ii = 1+1 = 2 ]
		ii = (1-i%2+i/2)*(P/3)+i*(1-(P/3)); // swap x<->y if P >= 3
		for(int j=0;j<3;j++) fM[3*ii+j] = M[ 3*((i+P)%3) + j]; // ii swaps xy, (i+P)%3 performs cyclic permutation P%3 times
	}
"""

def opMat(op=0):
	"""
	Reflections
	Diagonal elements 0,4,8
	R = [1-3] : M[(R-1)*4] = -1 														==>    M[4*((R-1+i)%3)] = -1 for i = 0
	R = [4-6] : M[(R-4)*4] = -1 and M[4*(R-4+1)%3] = -1    ==>    M[4*((R-1+i)%3)] = -1 for i = 0,1
	R = 7 : M[0] = -1 and M[4] = -1 and M[8] = -1					==>    M[4*((R-1+i)%3)] = -1 for i = 0,1,2
	Permutations
	i->ii     P<3: 012->012     else:   012->102   (swap x<->y if P >= 3)
	P<3:   P/3 =0, ii = i
	P>=3:  P/3 =1, ii = 1-i%2 + i/2  [  i<2:   ii = 1-i  (01->10)      i==2:   ii = 1+1 = 2 ]
	"""
	op = int(op)
	c = np.zeros(9)
	M = np.zeros((3,3))
	c[0], c[4], c[8] = 1.0, 1.0, 1.0
	R = op//6;
	P = op%6;
	ii = (R>0)*(1+(R-1)//3); # R->ii : 0->0, [1-3]->1, [4-6]->2, 7->3

	for i in range(ii):
		c[4*((R-1+i)%3)] = -1.;

	for i in range(3):
		ii = (1-i%2+i//2)*(P//3)+i*(1-(P//3)); # swap x<->y if P >= 3

		for j in range(3):
			M[ii][j] = c[ 3*((i+P)%3) + j]; # ii swaps xy, (i+P)%3 performs cyclic permutation P%3 times
	return M

def findMat(M):
	for s in range(48):
		d = np.linalg.norm(M-opMat(s))
		if d<1.0e-10:
			return s
	print("ERROR")

"""
table = np.zeros((48,48),int)
for s in range(48):
	A = opMat(s)
	for t in range(48):
		B = opMat(t)
		C = A.dot(B)
		table[s][t] = findMat(C)
np.savetxt("defdif/multtable",table,fmt="%d")
"""

def opMult(op1,op2):
	# M(op1) * M(op2) = M(opMult)
	return multtable[op1%48][op2%48]

def opInv(op):
	return (multtable[:,op]==0).argmax()
