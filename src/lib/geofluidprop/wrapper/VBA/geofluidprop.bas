' Example code for use in VBA
Attribute VB_Name = "geofluidprop"
#If Win64 Then
    Declare PtrSafe Function iapws95_get_rhoc Lib "geofluidprop.dll" () As Double
    Declare PtrSafe Function iapws92_sat_p_T Lib "geofluidprop.dll" (ByVal T As Double) As Double
#Else
    Declare Function iapws95_get_rhoc Lib "geofluidprop.dll" () As Double
    Declare Function iapws92_sat_p_T Lib "geofluidprop.dll" (ByVal T As Double) As Double
#End If

Public Function iapws92_pSatT(ByVal T As Double) As Double
    p = iapws92_sat_p_T(T)
    iapws92_pSatT = p * 0.000001
End Function
