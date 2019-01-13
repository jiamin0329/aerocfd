# The Strategy of Updating Buffer Block



Step 1: Load data from main block into buffer_send

Step 2: Send buffer_send to specific thread using mpi function

Step 3: Receive data which is sent to local thread using buffer_recv

Step 4: MPI_wait is implemented for MPI purpose

Step 5: Read data from buffer_recv

Step 6: Update buffer region



```mermaid
graph TD
	subgraph Update Buffer Region
		LoadDataToBufferSend-->SendBufferDataToSpecificThread
		SendBufferDataToSpecificThread-->ReceiveDataAtSpecificTHread
		ReceiveDataAtSpecificTHread-->MPI_wait
		MPI_wait-->ReadDataFromBufferReceive
		ReadDataFromBufferReceive-->UpdateBuffer
	end 
```

  

```mermaid
graph LR
	VariablesToTransfer(Data to Transfer)-->GeomInfo
	subgraph Geometry Info
		GeomInfo-->Spacing(Spacing, SpaceingI/J/K)
		GeomInfo-->Distance(WallDistance)
		GeomInfo-->inv_j(inv_jac)
		GeomInfo-->DxiDx(dxidx)
		DxiDx-->dxidx(dxidx,dxidy,dxidz)
   	 	DxiDx-->detadx(detadx,detady,detadz)
    	DxiDx-->dzetadx(dzetadx,dzetady,dzetadz)
    end
    VariablesToTransfer-->NSvar(NS Variables)
    subgraph NS Equation Variables
    	NSvar-->conVar(Conservative Variables)
    	NSvar-->tur(Variables in Turbulence Equation)
    end 
```



