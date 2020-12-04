----------------------------------------------------------------------------
-- Session library.
--
-- @release $Id: session.lua,v 1.29 2007/11/21 16:33:20 carregal Exp $
----------------------------------------------------------------------------

local lfs = require"lfs"
local serialize = require"cgilua.serialize"

local assert, error, ipairs, _G, loadfile, next, tostring, type = assert, error, ipairs, _G, loadfile, next, tostring, type
local format, gsub, strfind, strsub = string.format, string.gsub, string.find, string.sub
local tinsert = table.insert
local _open = io.open
local remove, time = os.remove, os.time
local mod, rand, randseed = math.mod, math.random, math.randomseed
local attributes, dir, mkdir = lfs.attributes, lfs.dir, lfs.mkdir

module ("cgilua.session")

local INVALID_SESSION_ID = "Invalid session identification"

----------------------------------------------------------------------------
-- Internal state variables.
local root_dir = nil
local timeout = 10 * 60 -- 10 minutes

--
-- Checks identifier's format.
--
local function check_id (id)
	return id and (strfind (id, "^%d+$") ~= nil)
end

--
-- Produces a file name based on a session.
-- @param id Session identification.
-- @return String with the session file name.
--
local function filename (id)
	return format ("%s/%s.lua", root_dir, id)
end

----------------------------------------------------------------------------
-- Deletes a session.
-- @param id Session identification.
----------------------------------------------------------------------------
function delete (id)
	if not check_id (id) then
		return nil, INVALID_SESSION_ID
	end
	remove (filename (id))
end

--
-- Searches for a file in the root_dir
--
local function find (file)
	local fh = _open (root_dir.."/"..file)
	if fh then
		fh:close ()
		return true
	else
		return false
	end
end

--
-- Creates a new identifier.
-- @return New identifier.
--
local function new_id ()
	return rand (999999999)
end

----------------------------------------------------------------------------
-- Creates a new session identifier.
-- @return Session identification.
----------------------------------------------------------------------------
function new ()
	local id = new_id ()
	if find (id..".lua") then
		repeat
			id = new_id (id)
		until not find (id..".lua")
		randseed (mod (id, 999999999))
	end
	return id
end

----------------------------------------------------------------------------
-- Changes the session identificator generator.
-- @param func Function.
----------------------------------------------------------------------------
function setidgenerator (func)
	if type (func) == "function" then
		new_id = func
	end
end

----------------------------------------------------------------------------
-- Loads data from a session.
-- @param id Session identification.
-- @return Table with session data or nil in case of error.
-- @return In case of error, also returns the error message.
----------------------------------------------------------------------------
function load (id)
	if not check_id (id) then
		return nil, INVALID_SESSION_ID
	end
	local f, err = loadfile (filename (id))
	if not f then
		return nil, err
	else
		return f()
	end
end

----------------------------------------------------------------------------
-- Saves data to a session.
-- @param id Session identification.
-- @param data Table with session data to be saved.
----------------------------------------------------------------------------
function save (id, data)
	if not check_id (id) then
		return nil, INVALID_SESSION_ID
	end
	local fh = assert (_open (filename (id), "w+"))
	fh:write "return "
	serialize (data, function (s) fh:write(s) end)
	fh:close()
end

----------------------------------------------------------------------------
-- Removes expired sessions.
----------------------------------------------------------------------------
function cleanup ()
	local rem = {}
	local now = time ()
	for file in dir (root_dir) do
		local attr = attributes(root_dir.."/"..file)
		if attr and attr.mode == 'file' then
			if attr.modification + timeout < now then
				tinsert (rem, file)
			end
		end
	end
	for _, file in ipairs (rem) do
		remove (root_dir.."/"..file)
	end
end

----------------------------------------------------------------------------
-- Changes the session timeout.
-- @param t Number of seconds to maintain a session.
----------------------------------------------------------------------------
function setsessiontimeout (t)
	if type (t) == "number" then
		timeout = t
	end
end

----------------------------------------------------------------------------
-- Changes the session directory.
-- @param path String with the new session directory.
----------------------------------------------------------------------------
function setsessiondir (path)
	path = gsub (path, "[/\\]$", "")
	-- Make sure the given path is a directory
	if not attributes (path, "mode") then
		assert (mkdir (path))
	end
	-- Make sure it can create a new file in the given directory
	local test_file = path.."/".._G.cgilua.tmpname()
	local fh, err = _open (test_file, "w")
	if not fh then
		error ("Could not open a file in session directory: "..
			tostring(err), 2)
	end
	fh:close ()
	remove (test_file)
	root_dir = path
end

----------------------------------------------------------------------------
local ID_NAME = "cgilua session identification"
local id = nil

----------------------------------------------------------------------------
-- Destroys the session.
----------------------------------------------------------------------------
function destroy ()
	data = {} -- removes data from session table to avoid recreation by `close'
	delete (id)
end

----------------------------------------------------------------------------
-- Open user session.
-- This function should be called before the script is executed.
----------------------------------------------------------------------------
function open ()
	-- Redefine cgilua.mkurlpath to manage the session identification
	local mkurlpath = _G.cgilua.mkurlpath
	function _G.cgilua.mkurlpath (script, data)
		if not data then
			data = {}
		end
		data[ID_NAME] = id
		return mkurlpath (script, data)
	end

	cleanup()

	id = _G.cgilua.QUERY[ID_NAME] or new()
	if id then
		_G.cgilua.QUERY[ID_NAME] = nil
		_G.cgilua.session.data = load (id) or {}
	end
end

----------------------------------------------------------------------------
-- Close user session.
-- This function should be called after the script is executed.
----------------------------------------------------------------------------
function close ()
	if next (_G.cgilua.session.data) then
		save (id, _G.cgilua.session.data)
		id = nil
	end
end

local already_enabled = false
----------------------------------------------------------------------------
-- Enables the use of sessions.
-- This function must be called by every script that needs sessions.
-- It just calls the `open' function and register the `close' function
-- to be called at the end of the execution.
----------------------------------------------------------------------------
function _G.cgilua.enablesession ()
	if already_enabled then -- avoid misuse when a script calls another one
		return
	else
		already_enabled = true
	end
	open ()
	_G.cgilua.addclosefunction (close)
end

