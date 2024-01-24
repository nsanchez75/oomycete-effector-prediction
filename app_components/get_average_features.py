import app_components.FEAT as FEAT

def get_average_features(sequence, df=0, protID=0):
	'''
	Method: Fills feature columns in a dataframe with net averages

	Input:

		- sequence: amino acid string
		- df: dataframe
		- protID: sequence identifier
	'''

	# get net cumulative sums
	#print("calculating for", sequence)
	length = min(len(sequence), 900)

	gravy = hydrophobicity = exposed = disorder = bulkiness = interface = 0.0
	for ind, aa in enumerate(sequence):
		if ind == length: break

		if aa.upper() in FEAT.INTERFACE_DIC:
			gravy += FEAT.GRAVY_DIC[aa.upper()]
			hydrophobicity += FEAT.HYDRO_DIC[aa.upper()]
			exposed += FEAT.EXPOSED_DIC[aa.upper()]
			disorder += FEAT.DISORDER_DIC[aa.upper()]
			bulkiness += FEAT.BULKY_DIC[aa.upper()]
			interface += FEAT.INTERFACE_DIC[aa.upper()]

	# return averages
	return [gravy           / length,
			hydrophobicity 	/ length,
			exposed 		/ length,
			disorder 		/ length,
			bulkiness 		/ length,
			interface 		/ length]