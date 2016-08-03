function varargout = PartVisualizerFigure(varargin)
% PARTVISUALIZER MATLAB code for PartVisualizer.fig
%      PARTVISUALIZER, by itself, creates a new PARTVISUALIZER or raises the existing
%      singleton*.
%
%      H = PARTVISUALIZER returns the handle to a new PARTVISUALIZER or the handle to
%      the existing singleton*.
%
%      PARTVISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARTVISUALIZER.M with the given input arguments.
%
%      PARTVISUALIZER('Property','Value',...) creates a new PARTVISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PartVisualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PartVisualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PartVisualizer

% Last Modified by GUIDE v2.5 27-Jul-2016 17:34:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PartVisualizerFigure_OpeningFcn, ...
                   'gui_OutputFcn',  @PartVisualizerFigure_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PartVisualizer is made visible.
function PartVisualizerFigure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PartVisualizer (see VARARGIN)

% Get name of the dataset.
datasetName = varargin{1};

% Choose default command line output for PartVisualizer
handles.output = hObject;

% Create a class to keep the data in.
handles.visualizerData = PartVisualizer(datasetName);

% Update GUI.
handles = handles.visualizerData.UpdateGUI(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PartVisualizer wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PartVisualizerFigure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
clear handles;


% --- Executes on slider movement.
function actSlider_Callback(hObject, eventdata, handles)
% hObject    handle to actSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newValue = handles.actSlider.Value;
handles = handles.visualizerData.ChangeThreshold(newValue, handles);


% --- Executes during object creation, after setting all properties.
function actSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to actSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in levelMenu.
function levelMenu_Callback(hObject, eventdata, handles)
% hObject    handle to levelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns levelMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from levelMenu
newValue = handles.levelMenu.Value;
handles.visualizerData = handles.visualizerData.FillPartInfo(newValue, 1);
handles = handles.visualizerData.UpdateGUI(handles);


% --- Executes during object creation, after setting all properties.
function levelMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to levelMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in partMenu.
function partMenu_Callback(hObject, eventdata, handles)
% hObject    handle to partMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns partMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from partMenu
newValue = handles.partMenu.Value;
handles.visualizerData = handles.visualizerData.FillPartInfo(handles.visualizerData.layerID, newValue);
handles = handles.visualizerData.UpdateGUI(handles);

% --- Executes during object creation, after setting all properties.
function partMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to partMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in upButton.
function upButton_Callback(hObject, eventdata, handles)
% hObject    handle to upButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.visualizerData.instanceOffset > 1
   handles.visualizerData.instanceOffset = handles.visualizerData.instanceOffset - 2; 
   handles = handles.visualizerData.VisualizeInstances(handles);
end


% --- Executes on button press in downButton.
function downButton_Callback(hObject, eventdata, handles)
% hObject    handle to downButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxOffset = max(1, (round(handles.visualizerData.numberOfInstances/2) * 2 - 5));
if handles.visualizerData.instanceOffset < maxOffset
   handles.visualizerData.instanceOffset = handles.visualizerData.instanceOffset + 2; 
   handles = handles.visualizerData.VisualizeInstances(handles);
end


% --- Executes on button press in beginButton.
function beginButton_Callback(hObject, eventdata, handles)
% hObject    handle to beginButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.visualizerData.instanceOffset ~= 1
   handles.visualizerData.instanceOffset = int32(1);
   handles = handles.visualizerData.VisualizeInstances(handles);
end


% --- Executes on button press in endButton.
function endButton_Callback(hObject, eventdata, handles)
% hObject    handle to endButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxOffset = max(1, (round(handles.visualizerData.numberOfInstances/2) * 2 - 5));
if handles.visualizerData.instanceOffset ~= maxOffset
   handles.visualizerData.instanceOffset = maxOffset; 
   handles = handles.visualizerData.VisualizeInstances(handles);
end
