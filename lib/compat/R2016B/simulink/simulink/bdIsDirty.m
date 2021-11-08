function isDirty = bdIsDirty(bdnames)
%bdIsDirty - Returns whether or not a loaded block diagram has unsaved changes
%
%   isDirty = bdIsDirty(bdname)
%
% bdname can be a character vector, a cell array of character vectors,
% or a double array.  All character vectors must be the names of loaded
% block diagrams.  All doubles must be the handles of loaded block diagrams.
% It is an error to supply an invalid handle, or the handle of anything other
% than a block diagram.  Do not specify a path or handle to a block or
% subsystem, or a block diagram that is not loaded.
%
% isDirty is a logical array with one entry for each block diagram.  The
% logical value is true if the block diagram has been modified in memory
% since it was loaded or last saved, and false if there are no unsaved
% changes.

%   Copyright 2016-2017 The MathWorks, Inc.

dirty = get_param(bdnames,'Dirty');
isDirty = strcmp(dirty,'on');
if ~ischar(bdnames)
	isDirty = reshape(isDirty,size(bdnames));
end